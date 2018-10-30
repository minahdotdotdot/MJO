# MJO Simulation
# Author: Minah Yang <lucia.yang@colorado.edu>
# Copyright 2018. Not for public distribution.

struct MJO_params
    LL        :: Float64 #\hat{L}     ~ 10^6              (m)
    HH        :: Float64 #H           ~ 5000              (m)
    UU        :: Float64 #U           ~ 5                 (m/s)
    QQ        :: Float64 #\hat{Q}     ~ 0.05              (m)
    T_RC      :: Float64 #T_{RC}      ~ 16 days = 1382400 (s)
    T_Q       :: Float64 #\tau_q      ~ 96 hours = 345600 (s)

    g         :: Float64 # gravity on earth : 9.80665     (m/s^2)
    RE        :: Float64 # radius of earth  : 6,371,000   (m)

    AA        :: Float64 # \alpha \approx 1.0183548918336027
    BB        :: Float64 # \beta \approx 750.
    DD        :: Float64 # CM \delta \approx 1.1
    Qs        :: Float64 # moisture saturation 56.2 to 70. % by ND constant, Qs = 0.05
    B         :: Float64 # CM "b" \approx 11.4

    deg       :: Float64 # Size of grid in degrees (0.25)
    lat_range :: Array{Float64, 1} # start and end latitudes
    lon_range :: Array{Float64, 1} # start and end longitudes

    lat       :: Array{Float64, 1} # Latitude
    lon       :: Array{Float64, 1} # Longitude  
    y         :: Array{Float64, 1} # y
    delt_x    :: Float64 # 1/delt_x (to eliminate divisions)
    delt_y    :: Float64 # 1/delt_y (to eliminate divisions)

    Ro        :: Float64 # \frac{1}{Ro*y} = 4*pi*L^2/(U*day in seconds * RE)
    Fr        :: Float64 # 1/Fr^2 = gH/U^2
    BQH       :: Float64 # \frac{\hat{Q}}{H}
    Tratio    :: Float64 # :L / (U*T_RC)

    PP        :: Float64 # Normalizing constant for C&M precipitation function.
    KK        :: Float64 # KK*Δ(quantity), need to meet CFL : KK < delt_x^2/h_time

    MJO_params(
        LL, HH, UU, QQ, T_RC, T_Q,
        g, RE, 
        AA, BB, DD, Qs, B, 
        deg, lat_range, lon_range, 
        PP, h_time
        ) =
    new(copy(LL), copy(HH), copy(UU), copy(QQ), copy(T_RC), copy(T_Q),
        copy(g), copy(RE),
        copy(AA), copy(BB), copy(DD), copy(Qs)/QQ, copy(B),
        copy(deg), copy(lat_range), copy(lon_range),
        lat_range[1]-deg/2: deg: lat_range[2]+deg/2, lon_range[1]: deg: lon_range[2]-deg,
        pi*copy(RE)*(lat_range[1]-deg/2: deg: lat_range[2]+deg/2)/(180.0*copy(LL)),
        180.0*copy(LL)/(copy(deg)*pi*copy(RE)), 180.0*copy(LL)/(copy(deg)*pi*copy(RE)),
        4.0*pi*LL^2/(3600.0*24.0*UU*RE), g*HH/UU^2, BB*QQ/HH, LL/(UU*T_RC), 
        copy(PP), 0.005*(copy(deg)*pi*copy(RE))^2/(180.0*copy(LL))^2/h_time
        )
end

struct MJO_State
    m1     :: Array{Float64, 2}
    n1     :: Array{Float64, 2}
    m2     :: Array{Float64, 2}
    n2     :: Array{Float64, 2}
    h1     :: Array{Float64, 2}
    h2     :: Array{Float64, 2}
    q      :: Array{Float64, 2}

    MJO_State(m1, n1, m2, n2, h1, h2, q) =
    new(copy(m1), copy(n1), copy(m2), copy(n2), copy(h1), copy(h2), copy(q))
end

import Base: +, -, *, /, maximum, minimum, getproperty

function+(A::MJO_State, B::MJO_State)
    return MJO_State(
        A.m1+B.m1, 
        A.n1+B.n1, 
        A.m2+B.m2, 
        A.n2+B.n2, 
        A.h1+B.h1,
        A.h2+B.h2,
        A.q+B.q)
end

function-(A::MJO_State, B::MJO_State)
    return MJO_State(
        A.m1-B.m1, 
        A.n1-B.n1, 
        A.m2-B.m2, 
        A.n2-B.n2, 
        A.h1-B.h1,
        A.h2-B.h2,
        A.q-B.q)
end

function*(A::MJO_State, c::T) where T<:Real
    return MJO_State(
        c*A.m1, 
        c*A.n1, 
        c*A.m2, 
        c*A.n2, 
        c*A.h1,
        c*A.h2,
        c*A.q)
end

function*(c::T, A::MJO_State) where T<:Real
    return MJO_State(
        c*A.m1, 
        c*A.n1, 
        c*A.m2, 
        c*A.n2, 
        c*A.h1,
        c*A.h2,
        c*A.q)
end

function/(A::MJO_State, c::T) where T<:Real
    d = 1/c
    return MJO_State(
        d*A.m1, 
        d*A.n1, 
        d*A.m2, 
        d*A.n2, 
        d*A.h1,
        d*A.h2,
        d*A.q)
end
function istherenan(A::Array{T,2}) where T<:Real
    if true in isnan.(A)
        return true
    end
    return false
end

function istherenan(A::MJO_State)
    for f in fieldnames(MJO_State)[5:7]
        if istherenan(getproperty(A, f))==true
            return true
        else
            return false
        end
    end
end

function isthereinf(A::Array{T,2}) where T<:Real
    if true in isinf.(A)
        return true
    end
    return false
end

function isthereinf(A::MJO_State)
    for f in fieldnames(MJO_State)[5:7]
        if isthereinf(getproperty(A, f))==true
            return true
        else
            return false
        end
    end
end

function getproperty(evol::Array{MJO_State,1}, field::Symbol)
    onefield = Array{Array{Float64,2},1}(undef,length(evol))
    for i = 1 : length(evol)
        onefield[i] = getproperty(evol[i],field)[2:end,:]
    end
    return onefield
end

function maximum(evolfield::Array{Array{Float64,2},1})
    return maximum(maximum.(evolfield))
end

function maximum(evol::Array{MJO_State,1}, field::Symbol)
    M = maximum(getproperty(evol[1], field))
    for i = 2 : length(evol)
        potentialmax = maximum(getproperty(evol[i], field)[2:end-1, :])
        if M < potentialmax
            M = potentialmax
        end
    end
    return M
end

function minimum(evolfield::Array{Array{Float64,2},1})
    return minimum(minimum.(evolfield))
end

function minimum(evol::Array{MJO_State,1}, field::Symbol)
    M = minimum(getproperty(evol[1], field))
    for i = 2 : length(evol)
        potentialmin = minimum(getproperty(evol[i], field)[2:end-1, :])
        if M > potentialmin
            M = potentialmin
        end
    end
    return M
end

function elemdiv(A::Array{Array{T,2},1}, B::Array{Array{T,2},1}) where T<:Real
    if length(A)-length(B)!=0
        error("evolfields must be same length.")
    end
    for i = 1 : length(A)
        A[i] = A[i]./B[i]
    end
    return A
end

function whereNaN(evol::Array{MJO_State,1}, field::Symbol)
    i=0
    for i = 1 : length(evol)
        m,n=size(getproperty(evol[1], field))
        for j = 1 : m
            for k = 1 : n
                if isnan(getproperty(evol[i], field)[j,k])== true
                    return i
                end
            end
        end
    end
    return length(A)+1
end

function evolfindNaN(evol::Array{MJO_State,1})
    for f in fieldnames(MJO_State)
       print(whereNaN(evol, f), "\t")
    end
end 

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

####################################################################################################
function gen_params(h_time::Float64)
    return MJO_params(10.0^6, # LL
                    5000.0,            # HH
                    5.0,               # UU
                    0.05,              # QQ
                    1382400.0,         # T_RC
                    345600,            # T_Q
                    9.80665,           #  g
                    6371000,           # RE
                    1.0184,            # AA
                    750.0,             # BB
                    1.1,               # DD
                    .058,              # Qs
                    11.4,              # B
                    0.25,              # degree
                    [-20.0, 20.0],     # lat_range
                    [0.0, 360.0],      # lon_range
                    17500.0,           # PP
                    h_time             # time-step length
                    )
end

h_time = 0.00005
params=gen_params(h_time);
grid_y = length(params.lat);
grid_x = length(params.lon);
