include("../codes/mjo_a.jl");
include("../codes/time_step.jl")

using PyPlot, Printf, JLD2, FileIO

function savecontourmaps(evol::Array{MJO_State,1}, str::String; draw::Symbol=:contourf)
    for f in fieldnames(MJO_State)
        evolfield = 1
        if f == :m1 || f ==:n1
            evolfield = elemdiv(getproperty(evol, f)[1:end-1,2:end-1], getproperty(evol, :h1)[1:end-1,2:end-1])
        elseif f == :m2 || f ==:n2
            evolfield = elemdiv(getproperty(evol, f)[1:end-1,2:end-1], getproperty(evol, :h2)[1:end-1,2:end-1])
        else
            evolfield = getproperty(evol, f)[1:end-1,2:end-1]
        end
        minval = minimum(evolfield);
        maxval = maximum(evolfield);
        for j = 1 : length(evolfield)
            fig, ax = subplots(); 
            fig[:set_size_inches](13,2); 
            ax[:set_aspect]("equal");
            fig[:colorbar](
                ax[draw](
                    params.lon[2:end-1], 
                    params.lat[1:end-1],  
                    evolfield[j],
                    vmin=minval,
                    vmax=maxval,
                    cmap="inferno"#"PuOr"
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
            evolfield = getproperty(state,f)[1:end-1, 2:end-1]./getproperty(state,:h1)[1:end-1, 2:end-1];
        elseif f ==:m2 || f ==:n2
            evolfield = getproperty(state,f)[1:end-1, 2:end-1]./getproperty(state,:h2)[1:end-1, 2:end-1];
        else
            evolfield = getproperty(state,f)[1:end-1, 2:end-1];
        end
        fig, ax = subplots(); 
        fig[:set_size_inches](13,2); 
        ax[:set_aspect]("equal");
        fig[:colorbar](
            ax[draw](
                params.lon[2:end-1], 
                params.lat[1:end-1],  
                evolfield,
                cmap="inferno"#"PuOr"
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
end

#function genInitSr(seed::Int)
#    random


# This is unit test E. 
x = range(0, stop=params.lon[end]*pi/180, length=grid_x);
ISE = MJO_State(
          zeros(grid_y, grid_x),                    # m1
          zeros(grid_y, grid_x),                    # n1
          zeros(grid_y, grid_x),                    # m2
          zeros(grid_y, grid_x),                    # n2
          2 .+repeat(sin.(x)', grid_y,1), # h1
          2 .+repeat(sin.(x)', grid_y,1), # h2
          repeat(range(10.0^(-16.0), stop=1.1*params.Qs, length=grid_y), 1, grid_x), # q
    );
#=

y = range(0, stop = 2*pi, length = grid_y);
ISDict = Dict("A" => MJO_State( ##PP = 15000
          zeros(grid_y, grid_x),                                       # m1
          zeros(grid_y, grid_x),                                       # n1
          zeros(grid_y, grid_x),                                       # m2
          zeros(grid_y, grid_x),                                       # n2
          ones(grid_y, grid_x),                                        # h1
          ones(grid_y, grid_x),                                   # h2
          repeat(range(10.0^(-16.0), stop=1.1*params.Qs, length=grid_y), 1, grid_x) # q
    ), 
    "B" => MJO_State(
          repeat(range(1,stop=2, length=grid_y), 1, grid_x), # m1
          repeat(range(3,stop=4, length=grid_x)', grid_y,1), # n1
          repeat(range(1,stop=2, length=grid_y), 1, grid_x), # m2
          repeat(range(3,stop=4, length=grid_x)',grid_y, 1), # n2
          ones(grid_y, grid_x),                     # h1
          ones(grid_y, grid_x),                     # h2
          repeat(range(10.0^(-16.0), stop=1.1*params.Qs, length=grid_y), 1, grid_x) # q
    ),
    "C" => MJO_State( #PP = 5000
          repeat(sin.(x)', grid_y,1), # m1
          zeros(grid_y, grid_x),      # n1
          repeat(sin.(x)', grid_y,1), # m2
          zeros(grid_y, grid_x),      # n2
          ones(grid_y, grid_x),       # h1
          ones(grid_y, grid_x),       # h2
          repeat(range(10.0^(-16.0), stop=1.1*params.Qs, length=grid_y), 1, grid_x) # q
    ), 
    "D" => MJO_State(
          zeros(grid_y, grid_x),      # m1
          repeat(sin.(y), 1, grid_x), # n1
          zeros(grid_y, grid_x),      # m2
          repeat(sin.(y), 1, grid_x), # n2
          ones(grid_y, grid_x),       # h1
          ones(grid_y, grid_x),       # h2
          repeat(range(10.0^(-16.0), stop=0.2*params.Qs, length=grid_y), 1, grid_x) # q
    ),
    "E" => MJO_State(
          zeros(grid_y, grid_x),                    # m1
          zeros(grid_y, grid_x),                    # n1
          zeros(grid_y, grid_x),                    # m2
          zeros(grid_y, grid_x),                    # n2
          2 .+repeat(sin.(x)', grid_y,1), # h1
          2 .+repeat(sin.(x)', grid_y,1), # h2
          repeat(range(10.0^(-16.0), stop=1.1*params.Qs, length=grid_y), 1, grid_x), # q
    ),
    "F" => MJO_State( #PP = 17500
          repeat(sin.(x)', grid_y,1), # m1
          repeat(sin.(y), 1, grid_x), # n1
          repeat(sin.(x)', grid_y,1), # m2ev
          repeat(sin.(y), 1, grid_x), # n2
          repeat(range(1,stop=2,length=grid_y), 1, grid_x), # h1
          repeat(range(2,stop=1,length=grid_y), 1, grid_x), # h2
          repeat(range(10.0^(-16.0), stop=1.1*params.Qs, length=grid_y), 1, grid_x) # q
    )
) =#


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




