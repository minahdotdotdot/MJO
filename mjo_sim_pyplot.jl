include("mjo_a.jl");
include("time_step.jl")
include("smooth_data.jl")

using PyPlot, Printf

function savecontourmaps(evol::Array{MJO_State,1}, str::String; 
    draw::Symbol=:contourf
    aspect::String=true)
    for f in fieldnames(MJO_State)
        evolfield = 1
        cm = "PuOr"
        if f == :m1 || f ==:n1
            evolfield = elemdiv(getproperty(evol, f), getproperty(evol, :h1))
        elseif f == :m2 || f ==:n2
            evolfield = elemdiv(getproperty(evol, f), getproperty(evol, :h2))
        else
            evolfield = getproperty(evol, f)
            cc = "BuGn"
        end
        minval = minimum(evolfield);
        maxval = maximum(evolfield);
        for j = 1 : length(evolfield)
            fig, ax = subplots(); 
            fig[:set_size_inches](13,2); 
            if aspect==true
                ax[:set_aspect]("equal");
            end
            fig[:colorbar](
                ax[draw](
                    params.lon, 
                    params.lat[2:end-1],  
                    evolfield[j],
                    vmin=minval,
                    vmax=maxval,
                    cmap=cc
                )
            );
            savefig(
                "../movies/"*string(f)*"/"*str*string(j),
                pad_inches=.10, 
                bbox_inches="tight"
            )
            close(fig)
        end
    end
end

@inline function savecontour(state::MJO_State, ii::String; draw::Symbol= :contourf)
    for f in fieldnames(MJO_State)
        evolfield = 1
        if f == :m1 || f ==:n1
            evolfield = getproperty(state,f)[2:end-1, :]./getproperty(state,:h1)[2:end-1, :];
        elseif f ==:m2 || f ==:n2
            evolfield = getproperty(state,f)[2:end-1, :]./getproperty(state,:h2)[2:end-1, :];
        else
            evolfield = getproperty(state,f)[2:end-1, :];
        end
        fig, ax = subplots(); 
        fig[:set_size_inches](13,2); 
        ax[:set_aspect]("equal");
        fig[:colorbar](
            ax[draw](
                params.lon, 
                params.lat[2:end-1],  
                evolfield,
                cmap="PuOr"
            )
        );
        savefig(
            "../movies/"*string(f)*"/"*ii,
            pad_inches=.10, 
            bbox_inches="tight"
        )
        close(fig)
    end
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
    savecontour(initial_state, str*"1")
    for i = 2 : N+1
        dxdt(params, state, tend);
        if istherenan(tend)==true || isthereinf(tend)==true
            error("We've got a NaN at "*string(i)*"!!! \n")
        end
        state = state + h * tend;
        if rem(i,every)==1
            savecontour(state, str*string(1+div(i,every)))
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
    savecontour(initial_state, str*"1")
    for i = 2 : N+1
        RK4_one(state, tend, params, h)
        if istherenan(tend)==true || isthereinf(tend)==true
            error("We've got a NaN at "*string(i)*"!!! \n")
        end
        state = state + h * tend;
        if rem(i,every)==1
            savecontour(state, str*string(1+div(i,every)))
        end
    end
    return state
end


function genInitSr(stencil::Array{T,2}=zeros(0,0)) where T<:Real
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
end



# Should have net north momentum = 0 (integral/sum) 
# & zero at boundaries
# check periodic boundary and make sure initial conditions are periodic in x
# zero momenta and constant heightfields
# random moisture field (random, or random then smooth it. )
# show velocity instead of momenta DONE
# see how long it takes for interesting moisture dynamics
    # interesting == breaking code
    # or dynamics 
    # 20 days ish
# github
# future: RK4 (implicit maybe)

# thoughts: make P be a function of concentration (instead of total water)
# i.e. Q <- Q/h1, Qs a different appropriate parameter 
# try getting rid of div terms system of 3 equations




