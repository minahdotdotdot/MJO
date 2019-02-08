# MJO Simulation
# Author: Minah Yang <lucia.yang@colorado.edu>
# Copyright 2018. Not for public distribution.
# structs and functions defined for structs

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
        copy(PP), 0.042 
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

struct MJO_State_im
    m1     :: Array{Complex{Float64}, 2}
    n1     :: Array{Complex{Float64}, 2}
    m2     :: Array{Complex{Float64}, 2}
    n2     :: Array{Complex{Float64}, 2}
    h1     :: Array{Complex{Float64}, 2}
    h2     :: Array{Complex{Float64}, 2}

    MJO_State_im(m1, n1, m2, n2, h1, h2) =
    new(copy(m1), copy(n1), copy(m2), copy(n2), copy(h1), copy(h2))
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
    for f in fieldnames(MJO_State)
        if istherenan(getproperty(A, f))==true
            return true
        else
            return false
        end
    end
end

function istherenan(A::MJO_State_im)
    for f in fieldnames(MJO_State_im)
        if istherenan(imag(getproperty(A, f)))==true || istherenan(real(getproperty(A, f)))==true
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
    for f in fieldnames(MJO_State)
        if isthereinf(getproperty(A, f))==true
            return true
        else
            return false
        end
    end
end

function isthereinf(A::MJO_State_im)
    for f in fieldnames(MJO_State_im)
        if isthereinf(imag(getproperty(A, f)))==true || isthereinf(real(getproperty(A, f)))==true
            return true
        else
            return false
        end
    end
end

function getproperty(evol::Array{MJO_State,1}, field::Symbol)
    onefield = Array{Array{Float64,2},1}(undef,length(evol))
    for i = 1 : length(evol)
        onefield[i] = getproperty(evol[i],field)[2:end-1,:]
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



include("smooth_data.jl")

function genInitSr(stencil::Array{T,2}=zeros(0,0); scheme::String="ex") where T<:Real
    if scheme =="ex"
        if stencil==zeros(0,0) # q is random field
            return MJO_State(
                zeros(grid_y, grid_x),        #m1
                zeros(grid_y, grid_x),        #n1
                zeros(grid_y, grid_x),        #m2
                zeros(grid_y, grid_x),        #m2
                ones(grid_y, grid_x),         #h1
                ones(grid_y, grid_x),         #h2
                rand(grid_y, grid_x)          #q
                )
        else              # q is random field smoothed
            return MJO_State(
                zeros(grid_y, grid_x),        #m1
                zeros(grid_y, grid_x),        #n1
                zeros(grid_y, grid_x),        #m2
                zeros(grid_y, grid_x),        #m2
                ones(grid_y, grid_x),         #h1
                ones(grid_y, grid_x),         #h2
                smoother(rand(grid_y,grid_x), stencil) #q
                )
       
        end
    elseif scheme == "imex"
        if stencil==zeros(0,0) # q is random field
            return MJO_State(
                zeros(grid_y, grid_x),        #m1
                zeros(grid_y, grid_x),        #n1
                zeros(grid_y, grid_x),        #m2
                zeros(grid_y, grid_x),        #m2
                zeros(grid_y, grid_x),        #eta1
                zeros(grid_y, grid_x),        #eta2
                rand(grid_y, grid_x)          #q
                )
        else              # q is random field smoothed
            return MJO_State(
                zeros(grid_y, grid_x),        #m1
                zeros(grid_y, grid_x),        #n1
                zeros(grid_y, grid_x),        #m2
                zeros(grid_y, grid_x),        #m2
                zeros(grid_y, grid_x),        #eta1
                zeros(grid_y, grid_x),        #eta2
                smoother(rand(grid_y,grid_x), stencil) #q
                )
        end
    elseif scheme =="im"
        grid_x2 = Int(grid_x/2+1)
            return MJO_State_im(
                zeros(Complex{Float64}, grid_y-1, grid_x2),  #m1_hat
                zeros(Complex{Float64}, grid_y-1, grid_x2),  #n1_hat
                zeros(Complex{Float64}, grid_y-1, grid_x2),  #m2_hat
                zeros(Complex{Float64}, grid_y-1, grid_x2),  #m2_hat
                zeros(Complex{Float64}, grid_y-1, grid_x2),  #eta1_hat
                zeros(Complex{Float64}, grid_y-1, grid_x2),  #eta2_hat
            )
    end
end
