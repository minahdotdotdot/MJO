# MJO Simulation
# Author: Minah Yang <lucia.yang@colorado.edu>
# Copyright 2018. Not for public distribution.

include("mjo_a.jl")
global h_time = 0.0009
params=gen_params(h_time);
global grid_y = length(params.lat);
global grid_x = length(params.lon);
IC = genInitSr(scheme="imex");
IChat = genInitSr(scheme="im");


###############################################################################
# DCT/DST & FFT and their inverses
using FFTW
#   - state: physical dim
#   - statehat: Coeffs in Cos/Sin/Fourier Domain

function dcsft(state::MJO_State, statehat::MJO_State_im,
    fields=fieldnames(MJO_State)[1:6], grid_x::Int=1440)
    for qq in fields
        if qq in fields[[1,3,5,6]]   # m1, m2, h1, h2 
            # DCT then RFFT
            getproperty(statehat,qq)[:,:] = FFTW.rfft(
                vcat(
                FFTW.r2r(getproperty(state,qq), FFTW.REDFT10, 1),
                zeros(1,grid_x)
                ), 2)
        else                        # n1, n2
            # DST then RFFT
            getproperty(statehat,qq)[:,:] = FFTW.rfft(
                vcat(
                    zeros(1,grid_x),
                    FFTW.r2r(getproperty(state,qq), FFTW.RODFT10, 1)
                    ), 2)
        end
    end
    return statehat
end

function idcsft(state::MJO_State, statehat::MJO_State_im,
    fields=fieldnames(MJO_State)[1:6], grid_x::Int=1440)
    for qq in fields
        if qq in fields[[1,3,5,6]]   # m1, m2, h1, h2 
            # iRFFT then iDCT
            getproperty(state,qq)[:,:] = .5/grid_y*FFTW.r2r(
                real(
                    FFTW.irfft(getproperty(statehat,qq), grid_x, 2)
                    )[1:end-1,:],
                FFTW.REDFT01, 1)
        else                        # n1, n2
            # iRFFT then iDST
            getproperty(state,qq)[:,:] = .5/grid_y*FFTW.r2r(
                real(
                    FFTW.irfft(getproperty(statehat,qq), grid_x, 2)
                    )[2:end,:],
                FFTW.RODFT01, 1)
        end
    end
    return state
end
#=
using LinearAlgebra
xx = pi/180*repeat(range(0,stop=5,length=grid_y),1,grid_x).*repeat(params.lon', grid_y,);
yy = pi/180*repeat(range(0,stop=5,length=grid_x)',grid_y).*repeat(params.lat, 1,grid_x);
OO = (sin.(xx)+cos.(xx)).*cos.(yy)
EE = (sin.(xx)+cos.(xx)).*sin.(yy)
IC = MJO_State(OO,EE,OO,OO,EE,OO,rand(grid_y,grid_x));
IChat = genInitSr(scheme="im");
origIC = deepcopy(IC);
dcsft(IC, IChat); idcsft(IC, IChat); diffIC = origIC-IC;
for qq in fieldnames(MJO_State)
    print(qq, ": ", norm(getproperty(diffIC, qq),2), "\n")
end
=#

###############################################################################
# EXPLICIT, NONLINEAR TERMS
# Params:
#     - params: system parameters.
#     - state: system state for whole spatially-extended system.
#     - exout: tendency of system is written here.
# h is really eta in the imex scheme.(h_i = 1 + \eta_i)


@inline function P(
        LL::T, UU::T, QQ::T, B::T, Qs::T, q::T, T_Q::T, PP::T) where T<:Real
    return (8.0*LL)/(QQ*86400000.0*UU*PP)*(exp(B*q/Qs)-1.0)
end

## Precipitation and Radiative Cooling together (RC never appears on its own.)
@inline function P_RC(
        LL::T, UU::T, QQ::T, B::T, Qs::T, q::T, T_Q::T, BQH::T, Tratio::T, hh2::T, hh1::T,
        PP::T) where T<:Real
    value = BQH*P(LL, UU, QQ, B, Qs, q, T_Q,PP)
    if (hh2 - hh1) > T(0)
        value = value - Tratio * (hh2-hh1)
    end
    return value
end

@inline function div_flux(
        m::Array{T, 2}, n::Array{T, 2}, h::Array{T, 2}, q::Array{T, 2}, 
        ii::Int64, iii::Int64, iiii::Int64, jj::Int64,
        delt_x::T, delt_y::T) where T<:Real
    return (
        .5*delt_x*(
        (m[jj,ii]+m[jj,iiii+1])/(2+h[jj,ii]+h[jj,iiii+1])*(q[jj,ii]+q[jj,iiii+1])
        - (m[jj,iii-1]+m[jj,ii])/(2+h[jj,iii-1]+h[jj,ii])*(q[jj,iii-1]+q[jj,ii])
        ) +
        .5*delt_y*(
            (n[jj,ii]+n[jj+1,ii])/(2+h[jj,ii]+h[jj+1,ii])*(q[jj,ii]+q[jj+1,ii])
            -(n[jj-1,ii]+n[jj,ii])/(2+h[jj-1,ii]+h[jj,ii])*(q[jj-1,ii]+q[jj,ii])
            )
        )
end 

@inline function h_sum(h1, h2, jj, ii)
    return h1[jj,ii] + h2[jj,ii]
end

@inline function h_sum_a(h1, h2, a, jj, ii)
    return h1[jj,ii] + a*h2[jj,ii]
end

function EXNL(params::MJO_params, state::MJO_State, exout::MJO_State)
    delt_x = params.delt_x        :: Float64
    delt_y = params.delt_y        :: Float64
    LL     = params.LL            :: Float64
    UU     = params.UU            :: Float64
    QQ     = params.QQ            :: Float64
    B      = params.B             :: Float64
    Qs     = params.Qs            :: Float64
    BQH    = params.BQH           :: Float64
    Tratio = params.Tratio        :: Float64
    AA     = params.AA            :: Float64
    DD     = params.DD            :: Float64
    T_Q    = params.T_Q           :: Float64
    PP     = params.PP            :: Float64
    KK     = params.KK            :: Float64

    # Ghost cells

    # Iterate over longitudinal direction.
    for ii = 1 : length(params.lon)
        # Account for zonal periodicity
        iii = ii  # left index: use for ii-1  
        iiii = ii # right index: use for ii+1  

        #Boundary conditions in lon: lon[1] = 0 = 360, lon[end] = long[1440] = 359.75
        if ii == 1
            iii = 1441 # to satisfy iii-1 = 1440
        elseif ii == 1440
            iiii = 0 #to satisfy iiii+1 = 1
        end

        # Ghost cells
        set_ghost_cells(state, ii)
        #println("column: ",ii)

        # Iterate over latitudinal direction.
        for jj = 2:length(params.y)-1
            #println("row:",jj)
            value_P_RC = P_RC(LL, UU, QQ, B, Qs, state.q[jj,ii], T_Q, BQH, Tratio, state.h2[jj,ii], state.h1[jj,ii],PP)

            ### MOMENTUM

            out.m1[jj,ii] = (
                - div_flux(state.m1, state.n1, state.h1, state.m1, ii, iii, iiii, jj, delt_x, delt_y)
                + params.Ro*params.y[jj]*state.n1[jj,ii]                             #=+1/Ro*n1=#
                - params.Fr*(
                    state.h1[jj,ii] * .5*delt_x*(h_sum(state.h1,state.h2,jj,iiii+1)-h_sum(state.h1,state.h2,jj,iii-1))
                    )                                                               #=h_1∂x(h_1+h_2)=#
                - state.m1[jj,ii]/(1+state.h1[jj,ii])*value_P_RC
                )
            #println("m1 done.")

            out.n1[jj,ii] = (
                - div_flux(state.m1, state.n1, state.h1, state.n1, ii, iii, iiii, jj, delt_x, delt_y)
                - params.Ro*params.y[jj]*state.m1[jj,ii]                             #=-1/Ro*m1=#
                - params.Fr*(
                    state.h1[jj,ii] * .5*delt_y*(h_sum(state.h1,state.h2,jj+1,ii)-h_sum(state.h1,state.h2,jj-1,ii))
                    )                                                               #=h_1∂y(h_1+h_2)=#
                - state.n1[jj,ii]/(1+state.h1[jj,ii])*value_P_RC
                ) 
            #println("n1 done.")

            out.m2[jj,ii] = (
                - div_flux(state.m2, state.n2, state.h2, state.m2, ii, iii, iiii, jj, delt_x, delt_y)
                + params.Ro*params.y[jj]*state.n2[jj,ii]                             #=+1/Ro*n2=#
                - params.Fr*(
                    state.h2[jj,ii] * .5*delt_x*(h_sum_a(state.h1,state.h2,AA,jj,iiii+1) - h_sum_a(state.h1,state.h2,AA,jj,iii-1))
                    )                                                               #=h_1∂x(h_1+\alpha* h_2)=#
                + state.m1[jj,ii]/(1+state.h1[jj,ii])*value_P_RC
                )
            #println("m2 done.")


            out.n2[jj,ii] = (
                - div_flux(state.m2, state.n2, state.h2, state.n2, ii, iii, iiii, jj, delt_x, delt_y)
                - params.Ro*params.y[jj]*state.m2[jj,ii]                             #=-1/Ro*m2=#
                - params.Fr*(
                    state.h2[jj,ii] * .5*delt_y*(h_sum_a(state.h1,state.h2,AA,jj+1,ii)-h_sum_a(state.h1,state.h2,AA,jj-1,ii))
                    )                                                               #=h_1∂y(h_1+\alpha*h_2)=#
                + state.n1[jj,ii]/(1+state.h1[jj,ii])*value_P_RC
        )
            #println("n2 done.")

            ### MASS

            out.h1[jj,ii] = - value_P_RC
            #println("h1 done.")

            out.h2[jj,ii] = + value_P_RC
            #println("h2 done.")

            ### MOISTURE
            out.q[jj,ii] = (
                - div_flux(state.m1, state.n1, state.h1, state.q, ii, iii, iiii, jj, delt_x, delt_y)
                +(-1.0+Qs./(DD*state.q[jj,ii]))*P(LL,UU,QQ,B,Qs,state.q[jj,ii], T_Q,PP) #=\hat{P}(Q)=#
                + KK*diffusion(state.q, ii, iii, iiii, jj, delt_x, delt_y)
                )
            #println("q done.")
        end
    end
    return exout
end

function imexstep(state::MJO_State, RHShat::MJO_State_im, out::MJO_State, params::MJO_params, h_time::Float64)
    # Calculate RHS. 
    expstate = deepcopy(state)
    exout(params, state, expstate);
    AA = params.AA
    Fr = params.Fr

    # Into Fourier/Cos/Sin Space
    dcsft(state + h_time * expstate, RHShat)
    a = h_time^2 * Fr*(kx.^2 + ky.^2);
    # b = -(1./a + ((1+AA) .+ (AA-1)*a ))

    # Implicit solve
    y4        = RHShat.m2 + im * kx./ky.*RHShat.n2
    y5        = RHShat.h2 - im * h_time * kx.*y4;
    y6        = (RHShat.h1 
        - im * h_time * kx.*RHShat.m1
        - h_time * ky * RHShat.n1
        + Fr/(h_time)*(1 .+a)./ky
        - (1 .+ 1./a).* y5
        )

    outhat    = deepcopy(RHShat);
    outhat.h2 = -1./(1./a + ((1+AA) .+ (AA-1)*a )) .* y6;
    outhat.n2 = 1/(h_time)* ky./(kx.^2 + ky.^2).*(y6-outhat.h2);
    outhat.m2 = y4 - im kx./ky.* outhat.n2;
    outhat.h1 = 1/Fr/h_time./ky .*(
        RHShat.n2 
        + AA*h_time*Fr*ky.*y6
        -y5
        );
    outhat.n1 = RHShat.n1 + h_time*Fr*ky.*(y6 + RHShat.n2);
    outhat.m1 = RHShat.m1 - im * h_time*Fr * kx.*(y6 + RHShat.n2);

    # Into Physical Space
    idcsft(out, outhat)
    return out
end




