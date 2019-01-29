include("smooth_data.jl")
include("mjo_a.jl")
include("mjo_ex.jl")
include("mjo_sim_pyplot.jl")
include("mjo_imex2.jl")

######################
## Inline Functions ##
######################
@inline function no0rem(x::Int, y::Int)
    r = rem(x,y)
    if r == 0
        return y
    else
        return r
    end
end

######################
## Explicit Schemes ##
######################

function f_euler(initial_state:: MJO_State, params::MJO_params, h::Float64, N::Int, every::Int)
    tend = deepcopy(initial_state)
    evol = Array{MJO_State,1}(undef,div(N, every)+1)
    evol[1] = initial_state
    state = deepcopy(initial_state)
    for i = 2 : N+1
        dxdt(params, state, tend);
        if istherenan(tend) == true || isthereinf(tend) ==true
            print(i)
            return evol[1:div(i, every)]
        end
        state = state + h * tend;
        if rem(i,every) == 1
            evol[1+div(i, every)] = state
        end
        
    end
    return evol
end

@inline function RK4_step(state::MJO_State, 
    tend::MJO_State, 
    tendlist=0, 
    i::Int, 
    params::MJO_params; 
    bb::Float64, h_time::Float64, scheme::Function=EXNL)
    yn = deepcopy(state);
    scheme2(params, state, tend);
    k = h * tend; #k1
    yn += 1/6*k;
    scheme2(params, state + .5*k, tend);
    k = h * tend; #k2
    yn += 1/3*k;
    scheme2(params, state + .5*k, tend);
    k = h * tend; #k3
    yn += 1/3*k;
    scheme2(params, state + k, tend);
    k = h * tend; #k4
    return yn + (1/6)*k, 0
end

function RK4(initial_state::MJO_State, params::MJO_params, h::Float64, N::Int, every::Int)
    evol = Array{MJO_State,1}(undef, div(N, every)+1)
    evol[1] = initial_state
    state = deepcopy(initial_state)
    tend = deepcopy(initial_state)
    for i = 2 : N+1
        state = RK4_step(state, tend, params, h)
        if rem(i,every)==1
            if istherenan(state)==true || isthereinf(state)==true
                print(i)
                return evol[1:div(i, every)]
            end
            evol[1+div(i, every)] = state
        end
    end
    return evol
end

@inline function ab1_step(state::MJO_State, 
    update::MJO_State, 
    tendlist::Array{MJO_State,1}, 
    i::Int, 
    params::MJO_params; 
    bb::Float64, h_time::Float64)
    #i>5 is the true index number.
    #This function will calculate the state at t_i. 
    #state is the state at t_{i-1}
    update =  state + h_time * tendlist[1] 
    tendlist[2] = EXNL(params, update, state; bb=bb, h_time=h_time)
    return update, tendlist
end

@inline function ab2_step(state::MJO_State, 
    update::MJO_State, 
    tendlist::Array{MJO_State,1}, 
    i::Int, 
    params::MJO_params; 
    bb::Float64, h_time::Float64,
    xind::Array{Array{Int,1},1}=[[2,1],[1,2]])
    #i>5 is the true index number.
    #This function will calculate the state at t_i. 
    #state is the state at t_{i-1}
    update =  state + h_time * (
        3/2 * tendlist[xind[no0rem(i,2)][1]] 
        -1/2 * tendlist[xind[no0rem(i,2)][2]] 
        )
    tendlist[no0rem(i,2)] = EXNL(params, update, state; bb=bb, h_time=h_time)
    return update, tendlist
end

@inline function ab3_step(state::MJO_State, 
    update::MJO_State, 
    tendlist::Array{MJO_State,1}, 
    i::Int, 
    params::MJO_params; 
    bb::Float64, h_time::Float64,
    xind::Array{Array{Int,1},1}=[[3,2,1],[1,3,2],[2,1,3]])
    #i>5 is the true index number.
    #This function will calculate the state at t_i. 
    #state is the state at t_{i-1}
    update =  state + h_time * (
        23/12 * tendlist[xind[no0rem(i,3)][1]] 
        -16/12 * tendlist[xind[no0rem(i,3)][2]] 
        +5/12 * tendlist[xind[no0rem(i,3)][3]]
        )
    tendlist[no0rem(i,3)] = EXNL(params, update, state; bb=bb, h_time=h_time)
    return update, tendlist
end

@inline function ab4_step(state::MJO_State, 
    update::MJO_State, 
    tendlist::Array{MJO_State,1}, 
    i::Int, 
    params::MJO_params; 
    bb::Float64, h_time::Float64,
    xind::Array{Array{Int,1},1}=[[4,3,2,1],[1,4,3,2],[2,1,4,3],[3,2,1,4]])
    #i>5 is the true index number.
    #This function will calculate the state at t_i. 
    #state is the state at t_{i-1}
    update =  state + h_time * (
        55/24 * tendlist[xind[no0rem(i,4)][1]] 
        -59/24 * tendlist[xind[no0rem(i,4)][2]] 
        +37/24 * tendlist[xind[no0rem(i,4)][3]] 
        -9/24 * tendlist[xind[no0rem(i,4)][4]] 
        )
    tendlist[no0rem(i,4)] = EXNL(params, update, state; bb=bb, h_time=h_time)
    return update, tendlist
end

######################
## Implicit Solve ####
######################
@inline function imsolve(exstate::MJO_State,
    RHShat::MJO_State_im, 
    outhat::MJO_State_im,
    params::MJO_params, h_time::Float64,
    kx::Array{Complex{Float64},2}, ky::Array{Float64,2}, 
    a::Array{Float64,2}, b::Array{Float64,2}, d::Array{Float64,2},
    f::Array{Float64,2}, g::Array{Float64,2})
    # Into Fourier/Cos/Sin Space
    dcsft(exstate, RHShat)

    # Implicit solve
    # Applying Lower^{-1}
    outhat.m1[:,:] = a .* RHShat.m1;
    outhat.n1[:,:] = a .* RHShat.n1;
    outhat.h1[:,:] = b .*(RHShat.h1 - (1/params.Fr) * (
        kx .*outhat.m1 + ky.*outhat.n1)); # b is actually 1/b
    outhat.m2[:,:] = a .* (RHShat.m2 - kx .* outhat.h1);
    outhat.n2[:,:] = a .* (RHShat.n2 + ky .* outhat.h1);
    outhat.h2[:,:] = (
        d .* (RHShat.h2 - 
        (1/params.Fr)* (kx .* outhat.m2 + ky .* outhat.n2))
        ); # d is actually a/d

    # Backward Substitution.
    # outhat.h2[:,:] = d .*outhat.h2;
    outhat.n2[:,:] = outhat.n2 .+ ky .* f.* outhat.h2;
    outhat.m2[:,:] = outhat.m2 -  kx .* f.* outhat.h2;
    outhat.h1[:,:] = outhat.h1 - g .* outhat.h2;
    outhat.n1[:,:] = outhat.n1 + ky.*a.*(outhat.h1 + outhat.h2);
    outhat.m1[:,:] = outhat.m1 - kx.*a.*(outhat.h1 + outhat.h2);

    # Into Physical Space
    idcsft(exstate, outhat)
    return exstate
end
#########################
## Save Time Evolution ##
#########################

#=exscheme can be: 
forward-euler (ab1_step), 
explicit RK4 (RK4_step), and 
Adams-Bashford with n steps (abn_step) =#


function imex(N::Int, every::Int, h_time::Float64; 
    bb::Float64=0.042, multistep::Bool=true, step::Int=3, exscheme::Function,
    msfunc::Array{Function,1}=[ab1_step, ab2_step, ab3_step, ab4_step])
    params = gen_params(h_time);
    IC     = genInitSr(scheme="imex");
    IChat  = genInitSr(scheme="im");
    state  = deepcopy(IC);     exstate = deepcopy(IC);
    RHShat = deepcopy(IChat);  outhat  = deepcopy(IChat);
    bb     = bb*h_time; #input bb should be the actual diffusion constant: K = bb/h_time.
    kx, ky, a, b, d, f, g = imex_init(params, h_time, bb);
    tendlist = 0; start = 2
    evol     = Array{MJO_State,1}(undef, div(N, every)+1);
    evol[1] = IC;
    if multistep==true
        tendlist = Array{MJO_State,1}(undef,step)
        for i = 1 : step-1
            tendlist[i] = EXNL(params, state, exstate, bb=bb, h_time=h_time)
            exstate, tendlist = msfunc[i](tendlist, i)
            exstate.q[:,:] = exstate.q + sqrt(h_time)*4.0e-7*tanh(3.0*exstate.q).*randn(size(exstate.q))
            state = imsolve(exstate, RHShat, outhat, params, h_time, kx, ky, a, b, d, f, g)
        end
        tendlist[step] = EXNL(params, state, exstate, bb=bb, h_time=h_time)
        start = step+1
        exscheme = msfunc[step]
    end
    for i = start : N+1
        exstate, tendlist = exscheme(state, exstate, tendlist, i, params, bb=bb, h_time=h_time)
        exstate.q[:,:] = exstate.q + sqrt(h_time)*4.0e-7*tanh(3.0*exstate.q).*randn(size(exstate.q))
        state = imsolve(exstate, RHShat, outhat, params, h_time,kx, ky, a, b, d, f, g)
        if rem(i, every) ==1
            evol[1+div(i,every)] = state
        end
    end
    return evol
end

###################
## Save Pictures ##
###################

function imex_print(N::Int, every::Int, h_time::Float64, name::String; 
    bb::Float64=0.042, multistep::Bool=true, step::Int=3, exscheme::Function,
    msfunc::Array{Function,1}=[ab1_step, ab2_step, ab3_step, ab4_step])
    params = gen_params(h_time);
    IC     = genInitSr(scheme="imex");
    IChat  = genInitSr(scheme="im");
    state  = deepcopy(IC);     exstate = deepcopy(IC);
    RHShat = deepcopy(IChat);  outhat  = deepcopy(IChat);
    bb     = bb*h_time; #input bb should be the actual diffusion constant: K = bb/h_time.
    kx, ky, a, b, d, f, g = imex_init(params, h_time, bb);
    tendlist = 0; start = 2; pad=ceil(Int,log10(N/every));
    saveimshow(state, name*string(1, pad=pad))
    if multistep==true
        tendlist = Array{MJO_State,1}(undef,step)
        for i = 1 : step-1
            tendlist[i] = EXNL(params, state, exstate, bb=bb, h_time=h_time)
            exstate, tendlist = msfunc[i](tendlist, i)
            exstate.q[:,:] = exstate.q + sqrt(h_time)*4.0e-7*tanh(3.0*exstate.q).*randn(size(exstate.q))
            state = imsolve(exstate, RHShat, outhat, params, h_time, kx, ky, a, b, d, f, g)
        end
        tendlist[step] = EXNL(params, state, exstate, bb=bb, h_time=h_time)
        start = step+1
        exscheme = msfunc[step]
    end
    for i = start : N+1
        exstate, tendlist = exscheme(state, exstate, tendlist, i, params, bb=bb, h_time=h_time)
        exstate.q[:,:] = exstate.q + sqrt(h_time)*4.0e-7*tanh(3.0*exstate.q).*randn(size(exstate.q))
        state = imsolve(exstate, RHShat, outhat, params, h_time,kx, ky, a, b, d, f, g)
        if rem(i, every) ==1
            saveimshow(state, name*string(1+div(i,every), pad=pad))
        end
    end
    return evol
end

