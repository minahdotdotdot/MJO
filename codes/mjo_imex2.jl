# MJO Simulation
# Author: Minah Yang <lucia.yang@colorado.edu>
# Copyright 2018. Not for public distribution.

include("mjo_a.jl")
global h_time = 0.0009 # = 3min
params=gen_params(h_time);
global grid_y = length(params.lat);
global grid_x = length(params.lon);
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
                    FFTW.r2r(getproperty(state,qq)[2:end-1,:], FFTW.REDFT10, 1),
                    zeros(1,grid_x)
                    ), 2)
        else                        # n1, n2
            # DST then RFFT
            getproperty(statehat,qq)[:,:] = FFTW.rfft(
                vcat(
                    zeros(1,grid_x),
                    FFTW.r2r(getproperty(state,qq)[2:end-1,:], FFTW.RODFT10, 1)#[1:end-1, :]
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
            getproperty(state,qq)[2:end-1,:] = .5/(grid_y-2)*FFTW.r2r(
                real(
                    FFTW.irfft(getproperty(statehat,qq), grid_x, 2)
                    )[1:end-1,:],
                FFTW.REDFT01, 1)
        else                        # n1, n2
            # iRFFT then iDST
            getproperty(state,qq)[2:end-1,:] = .5/(grid_y-2)*FFTW.r2r(
                real(
                    FFTW.irfft(
                        getproperty(statehat,qq), grid_x, 2)
                    )[2:end,:],
                FFTW.RODFT01, 1)
        end
    end
    return state
end

using LinearAlgebra
#=
xx = pi/180*repeat(range(0,stop=5,length=grid_y),1,grid_x).*repeat(params.lon', grid_y,);
yy = pi/180*repeat(range(0,stop=5,length=grid_x)',grid_y).*repeat(params.lat, 1,grid_x);
OO = (sin.(144000*xx)+cos.(2000*xx)).*cos.(160*yy)
EE = (sin.(144000*xx)+cos.(2000*xx)).*sin.(160*yy)
IC = MJO_State(OO,EE,OO,OO,EE,OO,rand(grid_y,grid_x));
IChat = genInitSr(scheme="im");
origIC = deepcopy(IC);
dcsft(IC, IChat); 
idcsft(IC, IChat); 
diffIC = origIC-IC;
for qq in fieldnames(MJO_State)
    print(qq, ": ", maximum(abs.(getproperty(diffIC, qq)[2:end-1,:])), "\n")
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

@inline function P(
        LL::T, UU::T, QQ::T, B::T, Qs::T, q::Array{T,2}, T_Q::T, PP::T) where T<:Real
    return (8.0*LL)/(QQ*86400000.0*UU*PP)*(exp.((B/Qs)*q) .- 1.0)
end

## Precipitation and Radiative Cooling together (RC never appears on its own.)
@inline function P_RC(
        LL::T, UU::T, QQ::T, B::T, Qs::T, q::T, T_Q::T, BQH::T, Tratio::T, hh2::T, hh1::T,
        PP::T) where T<:Real
    value = BQH*P(LL, UU, QQ, B, Qs, q, T_Q,PP)
    #H2 = 2 - H1 
    if (hh2 - hh1) > T(0)
        value = value - Tratio * (hh2-hh1)
    end
    return value
end

@inline function div_flux(
        m::Array{T, 2}, n::Array{T, 2}, h::Array{T, 2}, q::Array{T, 2}, 
        ii::Int64, iii::Int64, iiii::Int64, jj::Int64,
        delt_x::T, delt_y::T; H::T=1.0) where T<:Real
    return (
        .5*delt_x*(
        (m[jj,ii]+m[jj,iiii+1])/(2*H+h[jj,ii]+h[jj,iiii+1])*(q[jj,ii]+q[jj,iiii+1])
        - (m[jj,iii-1]+m[jj,ii])/(2*H+h[jj,iii-1]+h[jj,ii])*(q[jj,iii-1]+q[jj,ii])
        ) +
        .5*delt_y*(
            (n[jj,ii]+n[jj+1,ii])/(2*H+h[jj,ii]+h[jj+1,ii])*(q[jj,ii]+q[jj+1,ii])
            -(n[jj-1,ii]+n[jj,ii])/(2*H+h[jj-1,ii]+h[jj,ii])*(q[jj-1,ii]+q[jj,ii])
            )
        )
end 

@inline function diffusion(q::Array{T,2}, ii::Int64, iii::Int64, iiii::Int64, jj::Int64, 
    delt_x::T, delt_y::T) where T <:Real
    return (
        delt_x^2*(q[jj,iiii+1]-2*q[jj,ii]+q[jj,iii-1])
        +delt_y^2*(q[jj+1,ii]-2*q[jj,ii]+q[jj-1,ii])
        )
end

@inline function set_ghost_cells(state, ii)
    #= state.h1[1,ii]   = 2.0*state.h1[2,ii] - state.h1[3,ii]
    state.h2[1,ii]   = 2.0*state.h2[2,ii] - state.h2[3,ii]
    state.h1[end,ii] = 2.0*state.h1[end-1,ii] - state.h1[end-2,ii]
    state.h2[end,ii] = 2.0*state.h2[end-1,ii] - state.h2[end-2,ii]

    state.n1[1,ii]   = state.h1[1,ii] / state.h1[2,ii] * state.n1[2,ii]
    state.n2[1,ii]   = state.h2[1,ii] / state.h2[2,ii] * state.n2[2,ii]
    state.n1[end,ii] = state.h1[end,ii] / state.h1[end-1,ii] * state.n1[end-1,ii]
    state.n2[end,ii] = state.h2[end,ii] / state.h2[end-1,ii] * state.n2[end-1,ii] =#

    #=Ghost cells impose meridional boundary coniditons located at:
    j=1.5 = "-20 deg" & j=end-.5 = "20 deg" =#

    #= n_i(boundary)=0 =#
    state.n1[1,ii]   = -state.n1[2,ii]
    state.n1[end,ii] = -state.n1[end-1,ii]
    state.n2[1,ii]   = -state.n2[2,ii]
    state.n2[end,ii] = -state.n2[end-1,ii]

    #= ∂yh_i(boundary)=0 =#
    state.h1[1,ii]   = state.h1[2,ii]
    state.h1[end,ii] = state.h1[end-1,ii]
    state.h2[1,ii]   = state.h2[2,ii]
    state.h2[end,ii] = state.h2[end-1,ii]

    #= ∂ym_i(boundary)=0 =#
    state.m1[1,ii]   = state.m1[2,ii]
    state.m1[end,ii] = state.m1[end-1,ii]
    state.m2[1,ii]   = state.m2[2,ii]
    state.m2[end,ii] = state.m2[end-1,ii]

    #= ∂yq(boundary)=0 =#
    state.q[1,ii]   = state.q[2,ii]
    state.q[end,ii] = state.q[end-1,ii]
end

@inline function h_sum(h1, h2, jj, ii)
    return h1[jj,ii] + h2[jj,ii]
end

@inline function h_sum_a(h1, h2, a, jj, ii)
    return h1[jj,ii] + a*h2[jj,ii]
end

function EXNL(params::MJO_params, state::MJO_State, out::MJO_State; bb::Float64=0, h_time::Float64=0, H1::Float64=1.0)
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

    H2 = 2.0 - H1
    #Set diffusion constant for q equations. 
    KK = 0; 
    if bb==0 
        KK = params.KK            :: Float64
    else
        KK = bb/h_time
    end

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
            iiii = 0   # to satisfy iiii+1 = 1
        end

        # Ghost cells
        set_ghost_cells(state, ii)

        # Iterate over latitudinal direction.
        for jj = 2:length(params.y)-1
            value_P_RC = P_RC(LL, UU, QQ, B, Qs, state.q[jj,ii], T_Q, BQH, Tratio, state.h2[jj,ii], state.h1[jj,ii],PP)

            ### MOMENTUM

            out.m1[jj,ii] = ( -2.5*state.m1[jj,ii]
                - div_flux(state.m1, state.n1, state.h1, state.m1, ii, iii, iiii, jj, delt_x, delt_y, H=H1)
                + params.Ro*params.y[jj]*state.n1[jj,ii]                             #=+1/Ro*n1=#
                - params.Fr*(
                    state.h1[jj,ii] * .5*delt_x*(h_sum(state.h1,state.h2,jj,iiii+1)-h_sum(state.h1,state.h2,jj,iii-1))
                    )                                                               #=h_1∂x(h_1+h_2)=#
                - state.m1[jj,ii]/(H1+state.h1[jj,ii])*value_P_RC
                )

            out.n1[jj,ii] = ( -2.5*state.n1[jj,ii]
                - div_flux(state.m1, state.n1, state.h1, state.n1, ii, iii, iiii, jj, delt_x, delt_y, H=H1)
                - params.Ro*params.y[jj]*state.m1[jj,ii]                             #=-1/Ro*m1=#
                - params.Fr*(
                    state.h1[jj,ii] * .5*delt_y*(h_sum(state.h1,state.h2,jj+1,ii)-h_sum(state.h1,state.h2,jj-1,ii))
                    )                                                               #=h_1∂y(h_1+h_2)=#
                - state.n1[jj,ii]/(H1+state.h1[jj,ii])*value_P_RC
                )

            out.m2[jj,ii] = (
                - div_flux(state.m2, state.n2, state.h2, state.m2, ii, iii, iiii, jj, delt_x, delt_y, H=H2)
                + params.Ro*params.y[jj]*state.n2[jj,ii]                             #=+1/Ro*n2=#
                - params.Fr*(
                    state.h2[jj,ii] * .5*delt_x*(h_sum_a(state.h1,state.h2,AA,jj,iiii+1) - h_sum_a(state.h1,state.h2,AA,jj,iii-1))
                    )                                                               #=h_1∂x(h_1+\alpha* h_2)=#
                + state.m1[jj,ii]/(H1+state.h1[jj,ii])*value_P_RC
                )


            out.n2[jj,ii] = (
                - div_flux(state.m2, state.n2, state.h2, state.n2, ii, iii, iiii, jj, delt_x, delt_y, H=H2)
                - params.Ro*params.y[jj]*state.m2[jj,ii]                             #=-1/Ro*m2=#
                - params.Fr*(
                    state.h2[jj,ii] * .5*delt_y*(h_sum_a(state.h1,state.h2,AA,jj+1,ii)-h_sum_a(state.h1,state.h2,AA,jj-1,ii))
                    )                                                               #=h_1∂y(h_1+\alpha*h_2)=#
                + state.n1[jj,ii]/(H1+state.h1[jj,ii])*value_P_RC
                )

            ### MASS

            out.h1[jj,ii] = (- value_P_RC)

            out.h2[jj,ii] = (+ value_P_RC)

            ### MOISTURE
            out.q[jj,ii] = (
                - div_flux(state.m1, state.n1, state.h1, state.q, ii, iii, iiii, jj, delt_x, delt_y)
                +(-1.0+Qs./(DD*state.q[jj,ii]))*P(LL,UU,QQ,B,Qs,state.q[jj,ii], T_Q,PP) #=\hat{P}(Q)=#
                + KK*diffusion(state.q, ii, iii, iiii, jj, delt_x, delt_y)
                )
        end
    end
    return out
end

function feEXNL(params::MJO_params, state::MJO_State, tend::MJO_State, h_time::Float64; scheme2=true, bb::Float64=0)
    EXNL(params, state, tend, bb=bb, h_time=h_time);
    return state + h_time*tend
end

@inline function imex_init(params::MJO_params, h_time::Float64, bb::Float64; H1::Float64=1.0)
    grid_x2 = Int(grid_x/2+1);
    kx = ((params.LL/params.RE )* 
    repeat(range(0, stop=grid_x2-1)', grid_y-1,1));
    ky = ((9/2 * params.LL/params.RE)*
    repeat(range(0, stop=grid_y-2), 1,grid_x2));

    aa = H1 * (h_time^2 * params.Fr)*(kx.^2 + ky.^2);
    c = 1 .+ bb * (kx.^2 + ky.^2)
    ak = 1 ./ c; c = c.^2;
    b = 1 ./  (ak .*(aa + c)); # actually 1./b
    g = aa ./(aa + c)
    d = ak ./(1 .+ (-1 + params.AA)* ak.^2 .* aa .+ g) # actually ak ./d
    #f = (ak .*((-1 + params.AA)*aa .+ params.AA*c))./ (aa .+ c)
    f = ak .* (params.AA .- g)

    kx = H1 * (im*h_time*params.Fr) * kx;
    ky = H1 * (h_time*params.Fr) * ky
    return kx, ky, ak, b, d, f, g
end

@inline function genRandfield(;grid_y::Int64=162, grid_x::Int64=1440)
    R = randn(grid_y+4, grid_x); 
    R = 2*R[2:end-1,:] + R[1:end-2,:] + R[3:end,:]
    return 1/6*
    (   R[2:end-1,:]
        + hcat(R[2:end-1,2:end], R[2:end-1,1]) 
        + hcat(R[2:end-1, end], R[2:end-1, 1:end-1]) 
    )
end

@inline function AAval(;H1::Float64=1.0, g::Float64=9.80665, HH::Float64=5000.0)
    He = HH*H1*(1-.5*H1)
    return g*He/(g*He-22^2)
end

#=
IC = MJO_State(
rand(162,1440), rand(162,1440), 
rand(162,1440), rand(162,1440), 
rand(162,1440), rand(162,1440), 
rand(162,1440));
IChat = genInitSr(scheme="im")
=#




