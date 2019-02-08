#=
include("mjo_a.jl")
initial_state = MJO_State(
          zeros(grid_y, grid_x),                                       # m1
          zeros(grid_y, grid_x),                                       # n1
          zeros(grid_y, grid_x),                                       # m2
          zeros(grid_y, grid_x),                                       # n2
          ones(grid_y, grid_x),                                        # h1
          3.0* ones(grid_y, grid_x),                                   # h2
          repeat(range(10.0^(-16.0), stop=1.1*params.Qs, length=grid_y), 1, grid_x) # q
    );

h = 0.00005; N = 30 #Number of time-steps. 
global state = deepcopy(initial_state)
tend= deepcopy(state)
evol_m1= []; evol_n1=[]; evol_m2=[]; evol_n2=[]; evol_h1=[]; evol_h2=[]; evol_q=[];
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
=#
#=
using PyPlot
using PyCall
@pyimport matplotlib.animation as anim
#imshow(initial_state.q)
#fig = figure("heatmap")
#ax = PyPlot.axes()
#imshow(initial_state.q[161:-1:2,:])
#colorbar()
#, c=ColorGradient([:black, :white]))


#Construct Figure and Plot Data
fig = figure("moisture")
ax = PyPlot.axes()
global heatmap = ax[:imshow](zeros(size(initial_state.q)), aspect=8, origin="lower")

# Define the init function, which draws the first frame (empty, in this case)
function init()
    global heatmap
    heatmap[:set_data](Array{Float64,2}(undef, 162, 1442))
    return (heatmap, Union{})  # Union{} is the new word for None
end

# Animate draws the i-th frame, where i starts at i=0 as in Python.
function animate(i, var::Array{Array{Float64,2},1})
    global heatmap
    heatmap[:set_data](var[i])
    return (heatmap,Union{})
end

# Create the animation object by calling the Python function FuncAnimaton
myanim = anim.FuncAnimation(fig, animate, init_func=init, frames=100, interval=20)

# Convert it to an MP4 movie file and saved on disk in this format.
myanim[:save]("3Lines.mp4", bitrate=-1, extra_args=["-vcodec", "libx264", "-pix_fmt", "yuv420p"])
#myanim[:save]("test1.mp4", bitrate=-1, extra_args=["-vcodec", "libx264", "-pix_fmt", "yuv420p"])
    
# Function for creating an embedded video given a filename
function html_video(filename)
    open(filename) do f
        base64_video = base64encode(f)
        """<video controls src="data:video/x-m4v;base64,$base64_video">"""
    end
end

# Display the movie in a Julia cell as follows. Note it has animation controls for the user.
#display("text/html", html_video("3Lines.mp4"))
#display("text/html", html_video("test1.mp4"))
=#



