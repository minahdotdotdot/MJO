include("mjo_a.jl")
global h_time = 0.0009
params=gen_params(h_time);
global grid_y = length(params.lat);
 grid_x = length(params.lon);
IC = genInitSr();

# TOTALLY EXPLICIT SCHEME
##Craig and Mack Precipitation
#TODO(minah): add boolean/string to choose between {C&M ,Betts-Miller}?
#    Gregor: this would add undue overhead.
#            Instead, define a different function for B&M and pass the function as a parameter.
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
    #=return (
        .25*delt_x*(
            +(m[jj,ii]/h[jj,ii] + m[jj,iiii+1]/h[jj,iiii+1])*(q[jj,ii]+q[jj,iiii+1])
            - (m[jj,iii-1]/h[jj,iii-1] + m[jj,ii]/h[jj,ii])*(q[jj,iii-1]+q[jj,ii])
            ) + 
        .25*delt_y*(
            +(n[jj,ii]/h[jj,ii] + n[jj+1,ii]/h[jj+1,ii])*(q[jj,ii]+q[jj+1,ii])
            -(n[jj-1,ii]/h[jj-1,ii] + n[jj,ii]/h[jj,ii])*(q[jj-1,ii]+q[jj,ii])
            )
        )=#
    return (
        .5*delt_x*(
        (m[jj,ii]+m[jj,iiii+1])/(h[jj,ii]+h[jj,iiii+1])*(q[jj,ii]+q[jj,iiii+1])
        - (m[jj,iii-1]+m[jj,ii])/(h[jj,iii-1]+h[jj,ii])*(q[jj,iii-1]+q[jj,ii])
        ) +
        .5*delt_y*(
            (n[jj,ii]+n[jj+1,ii])/(h[jj,ii]+h[jj+1,ii])*(q[jj,ii]+q[jj+1,ii])
            -(n[jj-1,ii]+n[jj,ii])/(h[jj-1,ii]+h[jj,ii])*(q[jj-1,ii]+q[jj,ii])
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
    

@inline function h_sum(h1, h2, jj, ii)
    return h1[jj,ii] + h2[jj,ii]
end

@inline function h_sum_a(h1, h2, a, jj, ii)
    return h1[jj,ii] + a*h2[jj,ii]
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

# Compute tendency of the system.
# Params:
#     - params: system parameters.
#     - state: system state for whole spatially-extended system.
#     - out: tendency of system is written here.
function dxdt(params::MJO_params, state::MJO_State, out::MJO_State)
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
                - state.m1[jj,ii]/state.h1[jj,ii]*value_P_RC
                + KK*diffusion(state.m1, ii, iii, iiii, jj, delt_x, delt_y)
                )
            #println("m1 done.")

            out.n1[jj,ii] = (
                - div_flux(state.m1, state.n1, state.h1, state.n1, ii, iii, iiii, jj, delt_x, delt_y)
                - params.Ro*params.y[jj]*state.m1[jj,ii]                             #=-1/Ro*m1=#
                - params.Fr*(
                    state.h1[jj,ii] * .5*delt_y*(h_sum(state.h1,state.h2,jj+1,ii)-h_sum(state.h1,state.h2,jj-1,ii))
                    )                                                               #=h_1∂y(h_1+h_2)=#
                - state.n1[jj,ii]/state.h1[jj,ii]*value_P_RC
                + KK*diffusion(state.n1, ii, iii, iiii, jj, delt_x, delt_y)
                ) 
            #println("n1 done.")

            out.m2[jj,ii] = (
                - div_flux(state.m2, state.n2, state.h2, state.m2, ii, iii, iiii, jj, delt_x, delt_y)
                + params.Ro*params.y[jj]*state.n2[jj,ii]                             #=+1/Ro*n2=#
                - params.Fr*(
                    state.h2[jj,ii] * .5*delt_x*(h_sum_a(state.h1,state.h2,AA,jj,iiii+1) - h_sum_a(state.h1,state.h2,AA,jj,iii-1))
                    )                                                               #=h_1∂x(h_1+\alpha* h_2)=#
                + state.m1[jj,ii]/state.h1[jj,ii]*value_P_RC
                + KK*diffusion(state.m2, ii, iii, iiii, jj, delt_x, delt_y)
                )
            #println("m2 done.")


            out.n2[jj,ii] = (
                - div_flux(state.m2, state.n2, state.h2, state.n2, ii, iii, iiii, jj, delt_x, delt_y)
                - params.Ro*params.y[jj]*state.m2[jj,ii]                             #=-1/Ro*m2=#
                - params.Fr*(
                    state.h2[jj,ii] * .5*delt_y*(h_sum_a(state.h1,state.h2,AA,jj+1,ii)-h_sum_a(state.h1,state.h2,AA,jj-1,ii))
                    )                                                               #=h_1∂y(h_1+\alpha*h_2)=#
                + state.n1[jj,ii]/state.h1[jj,ii]*value_P_RC
                + KK*diffusion(state.n2, ii, iii, iiii, jj, delt_x, delt_y)
		)
            #println("n2 done.")

            ### MASS

            out.h1[jj,ii] = (
                -.5*delt_x * (state.m1[jj,iiii+1] - state.m1[jj,iii-1])             #=∂x m_1=#
                -.5*delt_y * (state.n1[jj+1,ii] - state.n1[jj-1,ii])                #=∂y n_1=#
                - value_P_RC
                + KK*diffusion(state.h1, ii, iii, iiii, jj, delt_x, delt_y)
                )
            #println("h1 done.")

            out.h2[jj,ii] = (
                -.5*delt_x * (state.m2[jj,iiii+1] - state.m2[jj,iii-1])             #=∂x m_2=#
                -.5*delt_y * (state.n2[jj+1,ii] - state.n2[jj-1,ii])                #=∂y n_2=#
                + value_P_RC
                + KK*diffusion(state.h2, ii, iii, iiii, jj, delt_x, delt_y)
                )
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
    return out
end