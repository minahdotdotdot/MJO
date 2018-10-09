function f_euler(initial_state:: MJO_State, params::MJO_params, h::Float64, N::Int, every::Int)
    tend = deepcopy(initial_state)
    evol = Array{MJO_State,1}(undef,div(N, every)+1)
    evol[1] = initial_state
    state = deepcopy(initial_state)
    for i = 2 : N+1
        dxdt(params, state, tend);
        if istherenan(tend) == true || isthereinf(tend) ==true
            return evol[1:div(i, every)]
        end
        state = state + h * tend;
        if rem(i,every) == 1
            evol[1+div(i, every)] = state
        end
        
    end
    return evol
end

## TODO: RK4? implicit schemes?