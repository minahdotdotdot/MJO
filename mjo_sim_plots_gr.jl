include("../codes/mjo_a.jl");
include("../codes/time_step.jl")
include("../codes/smooth_data.jl")

using Plots, Printf, JLD2, FileIO
gr()

function savecontourmaps(evol::Array{MJO_State,1}, str::String; fill_bool::Bool=true)
    for f in fieldnames(MJO_State)
        evolfield = 1
        if f == :m1 || f ==:n1
            evolfield = elemdiv(getproperty(evol, f)[2:end-1, :], getproperty(evol, :h1)[2:end-1, :])
        elseif f == :m2 || f ==:n2
            evolfield = elemdiv(getproperty(evol, f)[2:end-1, :], getproperty(evol, :h2)[2:end-1, :])
        else
            evolfield = getproperty(evol, f)[2:end-1, :]
        end
        minval = minimum(evolfield);
        maxval = maximum(evolfield);
        for j = 1 : length(evolfield)
            savefig(
                contour(
                    evolfield[j],
                    aspect_ratio=1,
                    clims=(minval, maxval),
                    fill=fill_bool
                    ),
                "../movies/"*string(f)*"/"*str*string(j)
                )
        end
    end
end

@inline function savecontour(state::MJO_State, ii::String; draw::Function=contour)
    for f in fieldnames(MJO_State)
        evolfield = 1
        if f == :m1 || f ==:n1
            evolfield = elemdiv(getproperty(evol, f)[2:end-1, :], getproperty(evol, :h1)[2:end-1, :])
        elseif f == :m2 || f ==:n2
            evolfield = elemdiv(getproperty(evol, f)[2:end-1, :], getproperty(evol, :h2)[2:end-1, :])
        else
            evolfield = getproperty(evol, f)[2:end-1, :]
        end
        GR.figure(size=(9,1))
        savefig(
            draw(
                evolfield,
                aspect_ratio=1,
                # fill=true
                ),
            "../movies/"*string(f)*"/"*ii
            )
    end
end

function f_euler_contour(
    initial_state:: MJO_State, 
    params::MJO_params, 
    h::Float64, 
    N::Int, 
    every::Int,
    str::String;
    draw::Function=contour
    )
    tend = deepcopy(initial_state)
    state = deepcopy(initial_state)
    savecontour(initial_state, str*"1")
    for i = 2 : N+1
        dxdt(params, state, tend);
        if istherenan(tend)==true || isthereinf(tend)==true
            return "We've got a NaN at "*string(i)*"!!! \n"
        end
        state = state + h * tend;
        if rem(i,every)==1
            savecontour(state, str*string(1+div(i,every)), draw=draw)
        end
    end
end

function genInitSr(smoothed::Bool=true)
    if smoothed==true 
        return MJO_State(
            zeros(grid_y, grid_x),        #m1
            zeros(grid_y, grid_x),        #n1
            zeros(grid_y, grid_x),        #m2
            zeros(grid_y, grid_x),        #m2
            ones(grid_y, grid_x),         #h1
            ones(grid_y, grid_x),         #h2
            smoother(rand(grid_y,grid_x)) #q
            )
    else
        return MJO_State(
            zeros(grid_y, grid_x),        #m1
            zeros(grid_y, grid_x),        #n1
            zeros(grid_y, grid_x),        #m2
            zeros(grid_y, grid_x),        #m2
            ones(grid_y, grid_x),         #h1
            ones(grid_y, grid_x),         #h2
            rand(grid_y, grid_x)          #q
            )
    end
end