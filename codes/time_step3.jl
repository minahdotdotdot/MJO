include("smooth_data.jl")
include("mjo_a.jl")
include("mjo_sim_pyplot.jl")
include("mjo_imex2.jl")
using Printf
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

##################################
## Single-step Explicit Schemes ##
##################################
# Adams-Bashforth 1 IS forward euler.
@inline function ab1_step(state::MJO_State, 
    update::MJO_State, 
    tendlist::Array{MJO_State,1}, 
    i::Int, 
    params::MJO_params; 
    bb::Float64, h_time::Float64, init::Bool=false, 
    H1::Float64=1.0,  msource::Array{Float64,1}=[0.21, 3.0])
    #i>5 is the true index number.
    #This function will calculate the state at t_i. 
    #state is the state at t_{i-1}
    update =  state + h_time * tendlist[1] 
    ind = 2
    if init==false
        ind = no0rem(i,1)
    end
    tendlist[ind] = EXNL(params, update, state; 
        bb=bb, h_time=h_time, H1=H1, msource=msource)
    return update, tendlist
end

@inline function RK4_step(state::MJO_State, 
    tend::MJO_State, 
    i::Int, 
    params::MJO_params; 
    bb::Float64, h_time::Float64, scheme2::Function=EXNL,
    H1::Float64=1.0, msource::Array{Float64,1}=[0.21, 3.0])
    yn = deepcopy(state);
    scheme2(params, state, tend, H1=H1, msource=msource);
    k = h * tend; #k1
    yn += 1/6*k;
    scheme2(params, state + .5*k, tend, H1=H1, msource=msource);
    k = h * tend; #k2
    yn += 1/3*k;
    scheme2(params, state + .5*k, tend, H1=H1, msource=msource);
    k = h * tend; #k3
    yn += 1/3*k;
    scheme2(params, state + k, tend, H1=H1, msource=msource);
    k = h * tend; #k4
    return yn + (1/6)*k, 0
end

@inline function ab2_step(state::MJO_State, 
    update::MJO_State, 
    tendlist::Array{MJO_State,1}, 
    i::Int, 
    params::MJO_params; 
    bb::Float64, h_time::Float64,
    xind::Array{Array{Int,1},1}=[[2,1],[1,2]],
    init::Bool=false,
    H1::Float64=1.0,  msource::Array{Float64,1}=[0.21, 3.0])
    #i>5 is the true index number.
    #This function will calculate the state at t_i. 
    #state is the state at t_{i-1}
    update =  state + h_time * (
        3/2 * tendlist[xind[no0rem(i,2)][1]] 
        -1/2 * tendlist[xind[no0rem(i,2)][2]] 
        )
    ind = 3
    if init==false
        ind = no0rem(i,2)
    end
    tendlist[ind] = EXNL(params, update, state; bb=bb, h_time=h_time, H1=H1, msource=msource)
    return update, tendlist
end

@inline function ab3_step(state::MJO_State, 
    update::MJO_State, 
    tendlist::Array{MJO_State,1}, 
    i::Int, 
    params::MJO_params; 
    bb::Float64, h_time::Float64,
    xind::Array{Array{Int,1},1}=[[3,2,1],[1,3,2],[2,1,3]],
    init::Bool=false,
    H1::Float64=1.0, msource::Array{Float64,1}=[0.21, 3.0])
    #i>5 is the true index number.
    #This function will calculate the state at t_i. 
    #state is the state at t_{i-1}
    update =  state + h_time * (
        23/12 * tendlist[xind[no0rem(i,3)][1]] 
        -16/12 * tendlist[xind[no0rem(i,3)][2]] 
        +5/12 * tendlist[xind[no0rem(i,3)][3]]
        )
    ind = 4;
    if init ==false
        ind = no0rem(i,3)
    end
    tendlist[ind] = EXNL(params, update, state; bb=bb, h_time=h_time, H1=H1, msource=msource)
    return update, tendlist
end

@inline function ab4_step(state::MJO_State, 
    update::MJO_State, 
    tendlist::Array{MJO_State,1}, 
    i::Int, 
    params::MJO_params; 
    bb::Float64, h_time::Float64,
    xind::Array{Array{Int,1},1}=[[4,3,2,1],[1,4,3,2],[2,1,4,3],[3,2,1,4]], 
    H1::Float64=1.0,  msource::Array{Float64,1}=[0.21, 3.0])
    #i>5 is the true index number.
    #This function will calculate the state at t_i. 
    #state is the state at t_{i-1}
    update =  state + h_time * (
        55/24 * tendlist[xind[no0rem(i,4)][1]] 
        -59/24 * tendlist[xind[no0rem(i,4)][2]] 
        +37/24 * tendlist[xind[no0rem(i,4)][3]] 
        -9/24 * tendlist[xind[no0rem(i,4)][4]] 
        )
    tendlist[no0rem(i,4)] = EXNL(params, update, state; bb=bb, h_time=h_time, H1=H1, msource=msource)
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
    a::Array{Float64,2}, b::Array{Float64,2}, c::Array{Float64,2}, d::Array{Float64,2},
    f::Array{Float64,2}, g::Array{Float64,2};
    H1::Float64=1.0)
    # Into Fourier/Cos/Sin Space
    dcsft(exstate, RHShat, params.grid_x)

    # Implicit solve
    # Applying Lower^{-1}
    #outhat.m1[:,:] = a .* RHShat.m1;
    #outhat.n1[:,:] = a .* RHShat.n1;
    outhat.h1[:,:] = b .*(RHShat.h1 - (1/params.Fr)* a.* (
        kx .* RHShat.m1 + ky .* RHShat.n1)); # b is actually 1/b
    outhat.m2[:,:] = c .* (RHShat.m2 - (2-H1)*kx .* outhat.h1); #c is actually 1/c
    outhat.n2[:,:] = c .* (RHShat.n2 + (2-H1)*ky .* outhat.h1);
    outhat.h2[:,:] = (
        d .* (RHShat.h2 - (1/params.Fr)*
            (kx .* outhat.m2 + ky .* outhat.n2))
        ); # d is actually 1/d

    # Backward Substitution.
    # outhat.h2[:,:] = d .*outhat.h2;
    outhat.n2[:,:] = outhat.n2 + (2-H1)*ky .* f.* outhat.h2;
    outhat.m2[:,:] = outhat.m2 - (2-H1)*kx .* f.* outhat.h2;
    outhat.h1[:,:] = outhat.h1 - g .* outhat.h2;
    outhat.n1[:,:] = a .*(RHShat.n1 + H1*ky .* (outhat.h1 + outhat.h2));
    outhat.m1[:,:] = a .*(RHShat.m1 - H1*kx .* (outhat.h1 + outhat.h2));

    # Into Physical Space
    idcsft(exstate, outhat, params.grid_x, params.grid_y)
    return exstate
end
#########################
## Save Time Evolution ##
#########################
function f_euler(initial_state:: MJO_State, params::MJO_params, h::Float64, N::Int, every::Int;
    H1::Float64=1.0,  msource::Array{Float64,1}=[0.21, 3.0])
    tend = deepcopy(initial_state)
    evol = Array{MJO_State,1}(undef,div(N, every)+1)
    evol[1] = initial_state
    state = deepcopy(initial_state)
    for i = 2 : N+1
        EXNL(params, state, tend, H1=H1, msource=msource);
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

function RK4(initial_state::MJO_State, params::MJO_params, h::Float64, N::Int, every::Int;
 H1::Float64=1.0, msource::Array{Float64,1}=[0.21, 3.0])
    evol = Array{MJO_State,1}(undef, div(N, every)+1)
    evol[1] = initial_state
    state = deepcopy(initial_state)
    tend = deepcopy(initial_state)
    for i = 2 : N+1
        state = RK4_step(state, tend, params, h, H1=H1)
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

#=exscheme can be: 
forward-euler (ab1_step), 
explicit RK4 (RK4_step), and 
Adams-Bashforth with n steps (abn_step) =#

function imex(N::Int, every::Int, h_time::Float64; 
    bb::Float64=0.005, multistep::Bool=true, step::Int=3, exscheme::Function=ab1_step,
    loc::String="../movies/",
    msfunc::Array{Function,1}=[ab1_step, ab2_step, ab3_step, ab4_step],
    params::MJO_params=gen_params(), IC=0, 
    NA::Float64=1.5, fr::Float64=0.5, msource::Array{Float64,1}=[0.21, 3.0])
    if IC==0
        IC = genInitSr(params, scheme="imex");
    end
    IChat  = genInitSr(params, scheme="im");
    state  = deepcopy(IC);     exstate = deepcopy(IC);
    RHShat = deepcopy(IChat);  outhat  = deepcopy(IChat);
    kx, ky, a, b, c, d, f, g = imex_init(params, h_time, bb, H1=H1, fr=fr);
    tendlist = Array{MJO_State,1}(undef, 1); start = 2;
    tendlist[1] = EXNL(params, state, exstate, bb=bb, h_time=h_time, H1=H1, msource=msource);
    if multistep==true
        tendlist = Array{MJO_State,1}(undef, step);
        tendlist[1] = EXNL(params, state, exstate, bb=bb, h_time=h_time, H1=H1, msource=msource);
        evol     = Array{MJO_State,1}(undef, div(N, every)+1);
        evol[1] = IC; tendlist[1] = EXNL(params, state, exstate, bb=bb, h_time=h_time, H1=H1, msource=msource);
        for i = 2 : step
            exstate, tendlist = msfunc[i-1](state, exstate, tendlist, i, params, bb=bb, h_time=h_time, init=true, H1=H1);
            #@printf("step %3d: maximum %4.2e \n",i, maximum(abs.(exstate.m1)))
            exstate.q[:,:] = exstate.q + sqrt(h_time)*NA*tanh.(3.0*exstate.q).*genRandfield(grid_y=params.grid_y, grid_x=params.grid_x) #7.45
            state = imsolve(exstate, RHShat, outhat, params, h_time, kx, ky, a, b, c, d, f, g, H1=H1)
            state.q[state.q .< 0.] .= 0.
        end
        start = step+1
    end
    for i = start : N+1
        exstate, tendlist = exscheme(state, exstate, tendlist, i, params, bb=bb, h_time=h_time, H1=H1, msource=msource)
        exstate.q[:,:] = exstate.q + sqrt(h_time)*NA*tanh.(3.0*exstate.q).*genRandfield(grid_y=params.grid_y, grid_x=params.grid_x) 
        #@printf("step %3d: maximum %4.2e \n",i, maximum(abs.(exstate.m1)))
        state = imsolve(exstate, RHShat, outhat, params, h_time,kx, ky, a, b, c, d, f, g, H1=H1)
        state.q[state.q .< 0.] .= 0.
        if rem(i, every) ==1
            if istherenan(state)==true || isthereinf(state)==true
                @printf("Nan alert!")
                return evol
            elseif minimum(state.h1) <= -1*H1
                @printf("H1 is too small.")
                return 0
            end
            evol[1+div(i,every)] = state
        end
    end
    return evol
end

###################
## Save Pictures ##
###################

function imex_print(N::Int, every::Int, h_time::Float64, name::String; 
    bb::Float64=0.005, multistep::Bool=true, step::Int=3, exscheme::Function=ab1_step,
    loc::String="../movies/",
    msfunc::Array{Function,1}=[ab1_step, ab2_step, ab3_step, ab4_step],
    params::MJO_params=gen_params(), IC=0, 
    NA::Float64=1.5, fr::Float64=0.5, msource::Array{Float64,1}=[0.21, 3.0],
    hov::Bool=false, everyH::Int)
    if IC==0
        IC = genInitSr(params, scheme="imex");
    end
    IChat  = genInitSr(params, scheme="im");
    state  = deepcopy(IC);     exstate = deepcopy(IC);
    RHShat = deepcopy(IChat);  outhat  = deepcopy(IChat);
    kx, ky, a, b, c, d, f, g = imex_init(params, h_time, bb, H1=H1, fr=fr);
    tendlist = Array{MJO_State,1}(undef, 1); start = 2;
    tendlist[1] = EXNL(params, state, exstate, bb=bb, h_time=h_time, H1=H1, msource=msource);
    pad=ceil(Int,log10(N/every));
    saveimshow(state, name*string(1, pad=pad),loc=loc, params=params)
    if hov==true
       newtxt!(state; name=name, loc=loc)
    end
    if multistep==true
        tendlist = Array{MJO_State,1}(undef, step);
        tendlist[1] = EXNL(params, state, exstate, bb=bb, h_time=h_time, H1=H1, msource=msource);
        for i = 2 : step
            exstate, tendlist = msfunc[i-1](state, exstate, tendlist, i, params, bb=bb, h_time=h_time, init=true, H1=H1);
            #@printf("step %3d: maximum %4.2e \n",i, maximum(abs.(exstate.m1)))
            exstate.q[:,:] = exstate.q + sqrt(h_time)*NA*tanh.(3.0*exstate.q).*genRandfield(grid_y=params.grid_y, grid_x=params.grid_x) #7.45
            state = imsolve(exstate, RHShat, outhat, params, h_time, kx, ky, a, b, c, d, f, g, H1=H1)
            state.q[state.q .< 0.] .= 0.
        end
        start = step+1
    end
    for i = start : N+1
        exstate, tendlist = exscheme(state, exstate, tendlist, i, params, bb=bb, h_time=h_time, H1=H1, msource=msource)
        exstate.q[:,:] = exstate.q + sqrt(h_time)*NA*tanh.(3.0*exstate.q).*genRandfield(grid_y=params.grid_y, grid_x=params.grid_x) 
        #@printf("step %3d: maximum %4.2e \n",i, maximum(abs.(exstate.m1)))
        state = imsolve(exstate, RHShat, outhat, params, h_time,kx, ky, a, b, c, d, f, g, H1=H1)
        state.q[state.q .< 0.] .= 0.
        if rem(i, every) ==1
            if istherenan(state)==true || isthereinf(state)==true
                @printf("Nan alert!")
                return 0
            elseif minimum(state.h1) <= -1*H1
                @printf("H1 is too small.")
                return 0
            elseif hov==true
                if rem(i, everyH)==1
                    addtxt!(state;name=name, loc=loc)
                end
            end
            saveimshow(state, name*string(1+div(i,every), pad=pad), loc=loc, params=params, H1=H1)
        end
    end
    newtxt!(state; name=name*"end", loc=loc, hov=false);#save last state.
end

using DelimitedFiles

@inline function newtxt!(state::MJO_State; name::String, loc::String="../movies/", hov::Bool=true)
    m =Int(size(state.q)[1]/2);
    for f in fieldnames(MJO_State)
        writedlm(loc*string(f)*"/"*name*".txt", getproperty(state,f)[m,:]')
        if hov==false
            writedlm(loc*string(f)*"/"*name*".txt", getproperty(state,f))
        end
        #(.5*(state[:f][m,:]+state[:f][m+1,:]))'
    end
end

@inline function addtxt!(state::MJO_State; name::String, loc::String="../movies/")
    m =Int(size(state.q)[1]/2);
    for f in fieldnames(MJO_State)
        open(loc*string(f)*"/"*name*".txt", "a") do io
            writedlm(io, getproperty(state, f)[m,:]')
        end
    end
end

function genInitSr(;name::String, loc::String="../movies/")
    return MJO_State(
        readdlm(loc*"m1/"*name*".txt"), #m1
        readdlm(loc*"n1/"*name*".txt"), #n1
        readdlm(loc*"m2/"*name*".txt"), #m2
        readdlm(loc*"n2/"*name*".txt"), #n2
        readdlm(loc*"h1/"*name*".txt"), #h1
        readdlm(loc*"h2/"*name*".txt"), #h2
        readdlm(loc*"q/"*name*".txt"),  #q
        )
end

function hovmollertxt(txtname::String, imagename::String; loc::String="../movies/", T=240)
    loc0 = deepcopy(loc)
    for f in fieldnames(MJO_State)
        cm = "bwr"; #"PuOr";
        if f == :q
            cm = "gist_ncar"#"BuGn";
        end
        if loc0 != "../movies/"
            loc = loc0*string(f)*"_"
        else
            loc = loc0*string(f)*"/"
        end
        fig, ax=subplots();
        getproperty(ax, :set_aspect)("equal");
        fig.colorbar(
            ax.imshow(
                readdlm(loc*txtname*".txt"),
                cmap=cm,
                extent=(params.lon_range[1], params.lon_range[2], T, 0)
            ),
            orientation="horizontal");
    savefig(
            loc*imagename*".png",
            pad_inches=.10,
            bbox_inches="tight"
        );
    close(fig)
    end
end