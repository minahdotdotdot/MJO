include("../codes/mjo_a.jl");
import Base.maximum, Base.minimum, Base.abs, Base./

function whereNaN(A::Array{Array{T,2},1}) where T<:AbstractFloat
	i=0
	for i = 1 : length(A)
		m,n=size(A[i])
		for j = 1 : m
			for k = 1 : n
				if isnan(A[i][j,k])== true
					return i
				end
			end
		end
	end
	return length(A)+1
end

function maximum(A::Array{Array{T,2},1}) where T<:AbstractFloat
	return maximum(maximum.(A))
end

function minimum(A::Array{Array{T,2},1})where T<:AbstractFloat
	return minimum(minimum.(A))
end

function abs(A::Array{Array{T,2},1}) where T<:AbstractFloat
	AA = deepcopy(A)
	for i = 1 : length(A)
		AA[i] = abs.(A[i])
	end
	return AA
end

function /(A::Array{Array{T,2},1}, c::TT) where {TT<:Real, T<:AbstractFloat}
	c = 1/c;
	for i = 1 : length(A)
		A[i] = c*A[i]
	end
	return A
end
#=
function saveheatmaps(evol, every::Int, str::String)
    N = length(evol);
    i = 1;
    while i*every <= N
        savefig(heatmap(evol[(i-1)*every+1], clim=(minimum(evol), 1)),"../movies/"*str*"_"*string(i))
        i+=1
    end
end
#, clim = [minimum(evol), maximum(evol)]

#include("../codes/unit_test.jl")
using Plots
pyplot()

initial_state = MJO_State(
          zeros(grid_y, grid_x),                                       # m1
          zeros(grid_y, grid_x),                                       # n1
          zeros(grid_y, grid_x),                                       # m2
          zeros(grid_y, grid_x),                                       # n2
          ones(grid_y, grid_x),                                        # h1
          3.0* ones(grid_y, grid_x),                                   # h2
          repeat(range(10.0^(-16.0), stop=1.1*params.Qs, length=grid_y), 1, grid_x) # q
);

h = 0.00005; N = 100 #Number of time-steps. 
global state = deepcopy(initial_state)
tend= deepcopy(state);
TT = Array{Float64,2}
evol_m1= Array{TT,1}(undef, 0); evol_n1=Array{TT,1}(undef, 0);  
evol_m2=Array{TT,1}(undef, 0); evol_n2=Array{TT,1}(undef, 0); 
evol_h1=Array{TT,1}(undef, 0); evol_h2=Array{TT,1}(undef, 0); 
evol_q=Array{TT,1}(undef, 0); 
for n = 1 : N
    dxdt(params, state, tend);
    global state = state + h * tend
    push!(evol_m1, state.m1)
    push!(evol_n1, state.n1)
    push!(evol_m2, state.m2)
    push!(evol_n2, state.n2)
    push!(evol_h1, state.h1)
    push!(evol_h2, state.h2)
    push!(evol_q, state.q)
end



#saveheatmaps(evol_h1, 2, "height1") 
=#
