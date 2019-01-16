include("smooth_data.jl")
include("mjo_a.jl")
include("mjo_ex.jl")
include("mjo_sim_pyplot.jl")
include("mjo_imex2.jl")


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

@inline function RK4_one(params::MJO_params, state::MJO_State, tend::MJO_State, h::Float64; 
    scheme2=dxdt::Function)
    #= # Option 1: Creates 4 variables, k1, k2, k3, k4, 
    # 3 scalar multiplications, 4 additions, 6 MJO_State variables at the end[2 initially +1 every stage]
    dxdt(params, state, tend);    
    k1 = h * tend;
    dxdt(params, state + .5*k1, tend);
    k = h * tend;
    dxdt(params, state + .5*k2, tend);
    k3 = h * tend;
    dxdt(params, state + k3, tend);
    k4 = h * tend;
    return state + (1/6)*(k1 + 2*k2 + 2*k3 + k4)=#

    # Option 2: Create one variable k(quarter storage compared to option 1), update yn at each stage
    # 4 scalar multiplications, 4 additions, 4 MJO_State varibles total [3 initially, +1: k updated.]
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
    return yn + (1/6)*k
end

function RK4(initial_state::MJO_State, params::MJO_params, h::Float64, N::Int, every::Int)
    evol = Array{MJO_State,1}(undef, div(N, every)+1)
    evol[1] = initial_state
    state = deepcopy(initial_state)
    tend = deepcopy(initial_state)
    for i = 2 : N+1
        state = RK4_one(state, tend, params, h)
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

@inline function no0rem(x::Int, y::Int)
    r = rem(x,y)
    if r == 0
        return y
    else
        return r
    end
end

function ab4_step(state::MJO_State, 
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

# scheme1: always gives the explicit update. 
#        : (NL^{(n+1)} = EXPLICITscheme(NL^{(n)})
#        : It is either feEXNL or RK4_one
# scheme2: necessary to tell RK4_one to use EXNL instead of dxdt. 
function imex_step(
    state::MJO_State, 
    exstate::MJO_State,
    RHShat::MJO_State_im, 
    outhat::MJO_State_im,
    params::MJO_params, h_time::Float64,
    kx::Array{Complex{Float64},2}, ky::Array{Float64,2}, 
    a::Array{Float64,2}, b::Array{Float64,2}, d::Array{Float64,2},
    f::Array{Float64,2}, g::Array{Float64,2};
    scheme1=feEXNL::Function, # Default is f_euler OR USE RK4_one
    scheme2=true,             # Default is f_euler OR USE EXNL
    bb::Float64=0.0
    )
    # Calculate RHS.
    exstate = scheme1(params, state, exstate, h_time; scheme2=scheme2, bb=bb);
    state = imsolve(exstate, RHShat, outhat, params, h_time,
    kx, ky, a, b, d, f, g)
    return state
end

function imex(
    IC::MJO_State, IChat::MJO_State_im, h_time::Float64, 
    N::Int, every::Int;
    scheme1=feEXNL::Function, #OR USE RK4_one
    scheme2=true, #OR USE EXNL
    bb::Float64=0
    )
    state           = deepcopy(IC);     exstate = deepcopy(IC);
    RHShat          = deepcopy(IChat);  outhat  = deepcopy(IChat);
    params = gen_params(h_time);
    bb = bb * (params.deg*pi*params.RE)^2/(180.0*params.LL)^2
    kx, ky, a, b, d, f, g = imex_init(params, h_time, bb);
    evol         = Array{MJO_State,1}(undef, div(N, every)+1)
    evol[1]      = IC;
    for i = 2 : N+1
        state = imex_step(
            state, exstate, RHShat, outhat, 
            params, h_time, kx, ky, a, b, d, f, g,
            scheme1=scheme1, scheme2=scheme2,
            bb=bb
            );
        if istherenan(state)==true||isthereinf(state)==true
            print(i)
            return evol[1:div(i,every)]
        end
        if rem(i, every) ==1
            #=if istherenan(state)==true||isthereinf(state)==true
                print(i)
                return evol[1:div(i,every)]
            end=#
            evol[1+div(i,every)] = state
            #end
        end
    end
    return evol
end

function imex_ab4(IC::MJO_State, IChat::MJO_State_im, h_time::Float64, 
    N::Int, every::Int;
    bb::Float64=0)
    state           = deepcopy(IC);     exstate = deepcopy(IC);
    RHShat          = deepcopy(IChat);  outhat  = deepcopy(IChat);

    params = gen_params(h_time);
    bb = bb * (params.deg*pi*params.RE)^2/(180.0*params.LL)^2
    kx, ky, a, b, d, f, g = imex_init(params, h_time, bb);

    #Do Adams-Bashford for s = 1, 2, 3. 
    evol        = Array{MJO_State,1}(undef, div(N, every)+1)
    tendlist    = Array{MJO_State,1}(undef,4)
    evol[1]     = IC;
    tendlist[1] = EXNL(params, state, exstate, bb=bb, h_time=h_time)
    #=2=# state = imsolve(state+h_time*tendlist[1], RHShat, outhat, params, h_time,
    kx, ky, a, b, d, f, g);
    tendlist[2] = EXNL(params, state, exstate, bb=bb, h_time=h_time)
    #=3=#state  = imsolve(state+ h_time*(3/2*tendlist[2]-1/2*tendlist[1]), 
        RHShat, outhat, params, h_time, kx, ky, a, b, d, f, g)
    tendlist[3] = EXNL(params, state, exstate, bb=bb, h_time=h_time)
    #=4=#state  = imsolve(state + h_time*(23/12*tendlist[3]-16/12*tendlist[2]+5/12*tendlist[1]), 
        RHShat, outhat, params, h_time, kx, ky, a, b, d, f, g)
    tendlist[4] = EXNL(params, state, exstate, bb=bb, h_time=h_time)
    # This code assumes: every>4. 
    for i = 5 : N+1
        exstate, tendlist = ab4_step(state, exstate, tendlist, i, params, bb=bb, h_time=h_time) #get exstate
        exstate.q[:,:] = exstate.q + sqrt(h_time)*4.0e-7*tanh(3.*exstate.q).*randn(size(exstate.q))
        state = imsolve(exstate, RHShat, outhat, params, h_time,kx, ky, a, b, d, f, g)
        if istherenan(state)==true||isthereinf(state)==true
            print(i)
            return evol[1:div(i,every)]
        end
        if rem(i, every) ==1
            evol[1+div(i,every)] = state
        end
    end
    return evol
end
function testimex_ab4(h_time::Float64, every::Int, name::String; bb::Float64=0)
    params = gen_params(h_time);
    IC    = genInitSr(scheme="imex");
    IChat = genInitSr(scheme="im");
    state           = deepcopy(IC);     exstate = deepcopy(IC);
    RHShat          = deepcopy(IChat);  outhat  = deepcopy(IChat);
    bb = bb * (params.deg*pi*params.RE)^2/(180.0*params.LL)^2 
    # (params.deg*pi*params.RE)^2/(180.0*params.LL)^2 /h_time is the CFL condition.
    # bb is some proportion of CFL condition s.t. bb/h_time == real diffusion constant

    kx, ky, a, b, d, f, g = imex_init(params, h_time, bb);
    @printf("i=   1: max = %4.2e, maxhat = %4.2e\n", 
                    maximum(abs.(state.m1)), maximum(norm.(outhat.m1)))
    N = Int(ceil(10*(365*24*60*60)/(h_time*2*10^5))); pad=ceil(Int,log10(N/every));
    #saveimshow(state, name*string(1, pad=pad))
    #Do Adams-Bashford for s = 1, 2, 3. 
    evol        = Array{MJO_State,1}(undef, div(N, every)+1)
    tendlist    = Array{MJO_State,1}(undef,4)
    evol[1]     = IC;
    tendlist[1] = EXNL(params, state, exstate, bb=bb, h_time=h_time)
    #=2=# state = imsolve(state+h_time*tendlist[1], RHShat, outhat, params, h_time,
    kx, ky, a, b, d, f, g);
    tendlist[2] = EXNL(params, state, exstate, bb=bb, h_time=h_time)
    #=3=#state  = imsolve(state+ h_time*(3/2*tendlist[2]-1/2*tendlist[1]), 
        RHShat, outhat, params, h_time, kx, ky, a, b, d, f, g)
    tendlist[3] = EXNL(params, state, exstate, bb=bb, h_time=h_time)
    #=4=#state  = imsolve(state + h_time*(23/12*tendlist[3]-16/12*tendlist[2]+5/12*tendlist[1]), 
        RHShat, outhat, params, h_time, kx, ky, a, b, d, f, g)
    tendlist[4] = EXNL(params, state, exstate, bb=bb, h_time=h_time)
    # This code assumes: every>4. 
    for i = 5 : N+1
        exstate, tendlist = ab4_step(state, exstate, tendlist, i, params, bb=bb, h_time=h_time) #get exstate
        exstate.q[:,:] = exstate.q + sqrt(h_time)*4.0e-7*tanh(3.*exstate.q).*randn(size(exstate.q))
        state = imsolve(exstate, RHShat, outhat, params, h_time, kx, ky, a, b, d, f, g)
        if rem(i, every) ==1
            if istherenan(state)==true || isthereinf(state)==true
                return i, state
            else
                @printf("i= %3d : max = %4.2e\n", 
                    i, maximum(abs.(state.m1)))
            end
            #saveimshow(state, name*string(1+div(i,every), pad=pad))
        end
    end
    return evol
end


function testimex_step(h_time::Float64, every::Int, name::String; bb::Float64=0)
    params = gen_params(h_time);
    IC    = genInitSr(scheme="imex");
    IChat = genInitSr(scheme="im");
    state           = deepcopy(IC);     exstate = deepcopy(IC);
    RHShat          = deepcopy(IChat);  outhat  = deepcopy(IChat);
    bb = bb * (params.deg*pi*params.RE)^2/(180.0*params.LL)^2 
    # (params.deg*pi*params.RE)^2/(180.0*params.LL)^2 /h_time is the CFL condition.
    # bb is some proportion of CFL condition s.t. bb/h_time == real diffusion constant

    kx, ky, a, b, d, f, g = imex_init(params, h_time, bb);
    @printf("i=   1: max = %4.2e, maxhat = %4.2e\n", 
                    maximum(abs.(state.m1)), maximum(norm.(outhat.m1)))
    N = Int(ceil(10*(365*24*60*60)/(h_time*2*10^5))); pad=ceil(Int,log10(N/every));
    saveimshow(state, name*string(1, pad=pad))
    for i = 2 : N # one year's time
        state = imex_step(
            state, exstate, RHShat, outhat, 
            params, h_time, kx, ky, a, b, d, f, g
            #, scheme1=RK4_one, scheme2=EXNL
            ,bb=bb
            );
        if rem(i, every) ==1
            if istherenan(state)==true || isthereinf(state)==true
                return i, state
            else
                @printf("i= %3d : max = %4.2e\n", 
                    i, maximum(abs.(state.m1)))
            end
            saveimshow(state, name*string(1+div(i,every), pad=pad))
        end
    end
    return outhat, state
end

bb=0.001
IC = genInitSr(scheme="imex"); IChat = genInitSr(scheme="im");
    state           = deepcopy(IC);     exstate = deepcopy(IC);
    RHShat          = deepcopy(IChat);  outhat  = deepcopy(IChat);
    params = gen_params(h_time);
    bb = bb * (params.deg*pi*params.RE)^2/(180.0*params.LL)^2
    kx, ky, a, b, d, f, g = imex_init(params, h_time, bb);

function testtime(state, exstate, RHShat, outhat, 
    params, h_time, kx, ky, a, b, d, f, g; bb)
    for i = 1 : 10
        outhat, state = imex_step(
                state, exstate, RHShat, outhat, 
                params, 0.0001, kx, ky, a, b, d, f, g, bb=bb
                );
    end
    return state
end

function f_euler_contour(
    initial_state:: MJO_State, 
    params::MJO_params, 
    h::Float64, 
    N::Int, 
    every::Int,
    str::String
    )
    tend = deepcopy(initial_state)
    state = deepcopy(initial_state)
    savecontourf(initial_state, str*"1")
    for i = 2 : N+1
        dxdt(params, state, tend);
        if istherenan(tend)==true || isthereinf(tend)==true
            error("We've got a NaN at "*string(i)*"!!! \n")
        end
        state = state + h * tend;
        if rem(i,every)==1
            savecontourf(state, str*string(1+div(i,every)))
        end
    end
    return state
end

function RK4_contour(
    initial_state:: MJO_State, 
    params::MJO_params, 
    h::Float64, 
    N::Int, 
    every::Int,
    str::String
    )
    tend = deepcopy(initial_state)
    state = deepcopy(initial_state)
    savecontourf(initial_state, str*"1")
    for i = 2 : N+1
        RK4_one(state, tend, params, h)
        if istherenan(tend)==true || isthereinf(tend)==true
            error("We've got a NaN at "*string(i)*"!!! \n")
        end
        state = state + h * tend;
        if rem(i,every)==1
            savecontourf(state, str*string(1+div(i,every)))
        end
    end
    return state
end

#=
bb = 10 .^(range(-4, stop=-1, length=4));
h_time = [0.0006, 0.0003, 0.0001, 0.00005, 0.00001];
for i = 1 : length(h_time)
    for j = 1 : length(bb)
    every = Int(floor(43200/(2*10^5*h_time[i])))
    @printf("================================================================\nh_time: %4.2e, bb: %3.2e, every: %d\n",h_time[i], bb[j], every);
    testimex_step(h_time[i], every, "halfday"*string(h_time[i]*2*10^5),bb= bb[j])
    end
end

bb = 10 .^(range(-4, stop=-1, length=4));
h_time = [0.0009, 0.0006, 0.0003, 0.0001, 0.00005, 0.00001];
for i = 1 : length(h_time)
    for j = 1 : length(bb)
        every = Int(floor(43200/(2*10^5*h_time[i])))
        print(every)
    end
end

=#

