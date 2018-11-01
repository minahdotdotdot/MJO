include("smooth_data.jl")
include("mjo_a.jl")
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

@inline function RK4_one(state::MJO_State, tend::MJO_State, params::MJO_params, h::Float64)
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
    dxdt(params, state, tend);
    k = h * tend; #k1
    yn += 1/6*k;
    dxdt(params, state + .5*k, tend);
    k = h * tend; #k2
    yn += 1/3*k;
    dxdt(params, state + .5*k, tend);
    k = h * tend; #k3
    yn += 1/3*k;
    dxdt(params, state + k, tend);
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


## TODO: RK4? implicit schemes?
