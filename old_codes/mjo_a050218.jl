# MJO Simulation
# Author: Minah Yang <lucia.yang@colorado.edu>
# Copyright 2018. Not for public distribution.

struct MJO_params
    LL      :: Float64 #\hat{L}     ~ 10^6              (m)
    HH      :: Float64 #H           ~ 5000              (m)
    UU      :: Float64 #U           ~ 5                 (m/s)
    QQ      :: Float64 #\hat{Q}     ~ 0.05              (m)
    T_RC    :: Float64 #T_{RC}      ~ 16 days = 1382400 (s)
    T_Q     :: Float64 #\tau_q      ~    96 hours = 345600 (s)

    g       :: Float64 # gravity on earth : 9.80665     (m/s^2)
    RE      :: Float64 # radius of earth  : 6,371,000   (m)

    AA      :: Float64 # \alpha \approx 1.0183548918336027
    BB      :: Float64 # \beta \approx 750.
    DD      :: Float64 # CM \delta \approx 1.1
    Qs      :: Float64 # moisture saturation 56.2 to 70. is reasonable dimensional input nondimensional output
    B       :: Float64 # CM "b" \approx 11.4

    Ro      :: Float64 # \frac{1}{Ro*y} = 4*pi*L^2/(U*day in seconds * RE)
    Fr      :: Float64 # 1/Fr^2 = gH/U^2
    BQH     :: Float64 # \frac{\hat{Q}}{H}
    Tratio  :: Float64 # :L / (U*T_RC)

    MJO_params(
        LL, HH, UU, QQ, T_RC, T_Q,
        g, RE, 
        AA, BB, DD, Qs, B, 
        grid_y, grid_x, delt_x, delt_y, lon, lat, y
        ) =
    new(copy(LL), copy(HH), copy(UU), copy(QQ), copy(T_RC), copy(T_Q),
        copy(g), copy(RE),
        copy(AA), copy(BB), copy(DD), copy(Qs)/QQ, copy(B),
        copy(grid_y), copy(grid_x), copy(delt_x), copy(delt_y), copy(lon), copy(lat), copy(y),
        4.*pi*LL^2/(3600.*24.*UU*RE), g*HH/UU^2, BB*QQ/HH, LL/(UU*T_RC))
end

struct MJO_State
    delt_x :: Float64 # 1/delt_x (to eliminate divisions)
    delt_y :: Float64 # 1/delt_y (to eliminate divisions)
    #TODO: we shouldn't use two copies of static things like the lat/long.
    lon    :: Array{Float64, 1}
    lat    :: Array{Float64, 1}

    y      :: Array{Float64, 1}

    m1     :: Array{Float64, 2}
    n1     :: Array{Float64, 2}
    m2     :: Array{Float64, 2}
    n2     :: Array{Float64, 2}
    h1     :: Array{Float64, 2}
    h2     :: Array{Float64, 2}
    q      :: Array{Float64, 2}

    MJO_State(delt_x, delt_y, lon, lat, y, m1, n1, m2, n2, h1, h2, q) =
    new(copy(delt_x), copy(delt_y), copy(lon), copy(lat), copy(y),
        copy(m1), copy(n1), copy(m2), copy(n2), copy(h1), copy(h2), copy(q))
end

##Craig and Mack Precipitation
#TODO(minah): add boolean/string to choose between {C&M ,Betts-Miller}?
#    Gregor: this would add undue overhead.
#            Instead, define a different function for B&M and pass the function as a parameter.
@inline function P{T<:Real}(
        LL::T, UU::T, QQ::T, B::T, Qs::T, q::T, T_Q::T)
    
    return (8.*LL)/(QQ*86400000.*UU)*(exp(B*q/Qs)-1.)
end

## Precipitation and Radiative Cooling together (RC never appears on its own.)
@inline function P_RC{T<:Real}(
        LL::T, UU::T, QQ::T, B::T, Qs::T, q::T, T_Q::T, BQH::T, Tratio::T, hh2::T, hh1::T)
    value = P(LL, UU, QQ, B, Qs, q, T_Q)
    if (hh2 - hh1) > T(0)
        value = value - Tratio * (hh2-hh1)
    end
    return value
end

@inline function div_flux{T<:Real}(
        m::Array{T, 2}, n::Array{T, 2}, h::Array{T, 2}, q::Array{T, 2}, ii::Int64, jj::Int64,
        delt_x::T, delt_y::T)
    #Account for zonal periodicity.
    iii = ii  #TODO(minah): annotate what this index is.
    iiii = ii  #TODO(minah): annotate what this index is.
    #TODO: is it possible to eliminate these branching points in this inner loop?
    if ii == 1
        iii = 1441 # to satisfy iii-1 = 1440
    elseif ii == 1440
        iiii = 0 #to satisfy iiii+1 = 1
    end

    #TODO: eliminate redundant division by delt
    return (
        .25*delt_x*(
            +(m[jj,ii]/h[jj,ii] + m[jj,iiii+1]/h[jj,iiii+1])*(q[jj,ii]+q[jj,iiii+1])
            - (m[jj,iii-1]/h[jj,iii-1] + m[jj,ii]/h[jj,ii])*(q[jj,iii-1]+q[jj,ii])
            ) + 
        .25*delt_y*(
            +(n[jj,ii]/h[jj,ii] + n[jj+1,ii]/h[jj+1,ii])*(q[jj,ii]+q[jj+1,ii])
            -(n[jj-1,ii]/h[jj-1,ii] + n[jj,ii]/h[jj,ii])*(q[jj-1,ii]+q[jj,ii])
            )
        )
end

@inline function h_sum(h1, h2, jj, ii)
    return h1[jj,ii] + h2[jj,ii]
end

@inline function h_sum_a(h1, h2, a, jj, ii)
    return h1[jj,ii] + a*h2[jj,ii]
end

@inline function set_ghost_cells(state, ii)
    state.h1[1,ii]   = 2.0*state.h1[2,ii] - state.h1[3,ii]
    state.h2[1,ii]   = 2.0*state.h2[2,ii] - state.h2[3,ii]
    state.h1[end,ii] = 2.0*state.h1[end-1,ii] - state.h2[end-2,ii]
    state.h2[end,ii] = 2.0*state.h2[end-1,ii] - state.h2[end-2,ii]

    state.n1[1,ii]   = state.h1[1,ii] / state.h1[2,ii] * state.n1[2,ii]
    state.n2[1,ii]   = state.h2[1,ii] / state.h2[2,ii] * state.n2[2,ii]
    state.n1[end,ii] = state.h1[end,ii] / state.h1[end-1,ii] * state.n1[end-1,ii]
    state.n2[end,ii] = state.h2[end,ii] / state.h2[end-1,ii] * state.n2[end-1,ii]
end

# Compute tendency of the system.
# Params:
#     - params: system parameters.
#     - state: system state for whole spatially-extended system.
#     - out: tendency of system is written here.
function dxdt(params::MJO_params, state::MJO_State, out::MJO_State)
    delt_x  =  state.delt_x         :: Float64
    delt_y  =  state.delt_y         :: Float64
    LL      =  params.LL            :: Float64
    UU      =  params.UU            :: Float64
    QQ      =  params.QQ            :: Float64
    B       =  params.B             :: Float64
    Qs      =  params.Qs            :: Float64
    BQH     =  params.BQH           :: Float64
    Tratio  =  params.Tratio        :: Float64
    AA      =  params.AA            :: Float64
    DD      =  params.DD            :: Float64
    T_Q     =  params.T_Q           :: Float64

    # Ghost cells

    # Iterate over longitudinal direction.
    for ii = 1 : length(state.lon)
        iii = ii  #Use for ii-1  #TODO(minah): annotate this index.
        iiii = ii #Use for ii+1  #TODO(minah): annotate this index.

        #Boundary conditions in lon: lon[1] = 0 = 360, lon[end] = long[1440] = 359.75
        if ii == 1
            iii = 1441 # to satisfy iii-1 = 1440
        elseif ii == 1440
            iiii = 0 #to satisfy iiii+1 = 1
        else
            set_ghost_cells(state, ii+1)
        end

        # Ghost cells
        set_ghost_cells(state,ii)

        #println("column: ",ii)

        # Iterate over latitudinal direction.
        for jj = 2:length(state.y)-1
            #println("row:",jj)
            value_P_RC = P_RC(LL, UU, QQ, B, Qs, state.q[jj,ii], T_Q, BQH, Tratio, state.h2[jj,ii], state.h1[jj,ii])

            ### MOMENTUM

            out.m1[jj,ii] = (
                + div_flux(state.m1, state.n1, state.h1, state.m1, ii, jj, delt_x, delt_y)
                + params.Ro*state.y[jj]*state.n1[jj,ii]                             #=1/Ro*n1=#
                - params.Fr*(
                    state.h1[jj,ii] * .5*delt_x*(h_sum(state.h1,state.h2,jj,iiii+1)-h_sum(state.h1,state.h2,jj,iii-1))
                    )                                                               #=h_1∂x(h_1+h_2)=#
                - state.m1[jj,ii]/state.h1[jj,ii]*value_P_RC
                )
            #println("m1 done.")

            out.n1[jj,ii] = (
                + div_flux(state.m1, state.n1, state.h1, state.n1, ii, jj, delt_x, delt_y)
                - params.Ro*state.y[jj]*state.m1[jj,ii]                             #=-1/Ro*m1=#
                - params.Fr*(
                    state.h1[jj,ii] * .5*delt_y*(h_sum(state.h1,state.h2,jj+1,ii)-h_sum(state.h1,state.h2,jj-1,ii))
                    )                                                               #=h_1∂y(h_1+h_2)=#
                - state.n1[jj,ii]/state.h1[jj,ii]*value_P_RC
                )
            #println("n1 done.")

            out.m2[jj,ii] = (
                + div_flux(state.m2, state.n2, state.h2, state.m2, ii, jj, delt_x, delt_y)
                + params.Ro*state.y[jj]*state.n2[jj,ii]                             #=1/Ro*n2=#
                - params.Fr*(
                    state.h1[jj,ii] * .5*delt_x*(h_sum_a(state.h1,state.h2,AA,jj,iiii+1) - h_sum_a(state.h1,state.h2,AA,jj,iii-1))
                    )                                                               #=h_1∂x(h_1+\alpha* h_2)=#
                + state.m1[jj,ii]/state.h1[jj,ii]*value_P_RC
                )
            #println("m2 done.")


            out.n2[jj,ii] = (
                + div_flux(state.m2, state.n2, state.h2, state.n2, ii, jj, delt_x, delt_y)
                - params.Ro*state.y[jj]*state.m2[jj,ii]                             #=-1/Ro*m2=#
                - params.Fr*(
                    state.h1[jj,ii] * .5*delt_y*(h_sum_a(state.h1,state.h2,AA,jj+1,ii)-h_sum_a(state.h1,state.h2,AA,jj-1,ii))
                    )                                                               #=h_1∂y(h_1+\alpha*h_2)=#
                + state.n1[jj,ii]/state.h1[jj,ii]*value_P_RC
                )
            #println("n2 done.")

            ### MASS

            out.h1[jj,ii] = (
                -.5*delt_x * (state.m1[jj,iiii+1] - state.m1[jj,iii-1])             #=∂x m_1=#
                -.5*delt_y * (state.n1[jj+1,ii] - state.n1[jj-1,ii])                #=∂y n_1=#
                - value_P_RC
                )
            #println("h1 done.")

            out.h2[jj,ii] = (
                -.5*delt_x * (state.m2[jj,iiii+1] - state.m2[jj,iii-1])             #=∂x m_2=#
                -.5*delt_y * (state.n2[jj+1,ii] - state.n2[jj-1,ii])                #=∂y n_2=#
                + value_P_RC
                )
            #println("h2 done.")

            ### MOISTURE
            out.q[jj,ii] = (
                + div_flux(state.m1, state.n1, state.h1, state.q, ii, jj, delt_x, delt_y)
                +(-1.+Qs./(DD*state.q[jj,ii]))*P(LL,UU,QQ,B,Qs,state.q[jj,ii], T_Q) #=\hat{P}(Q)=#
                )
            #println("q done.")
        end
    end
    return out
end

####################################################################################################

grid_y = 160  # Number of Rows (latitude) (jj)
grid_x = 1440 # Number of columns (longitude) (ii)

params = MJO_params(10.^6,    #LL
                    5000.,    #HH
                    5.,       #UU
                    0.05,     #QQ
                    1382400., #T_RC
                    345600,   #T_Q
                    9.80665,  # g
                    6371000,  #RE
                    1.0184,   #AA
                    750.,     #BB
                    1.1,      #DD
                    .058,     #Qs
                    11.4      #B
);







