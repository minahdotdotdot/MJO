# MJO
## mjo_a.jl
Sets up 2 structs:
 - MJO_State: fields are the 7 variables getfields(MJO_State) = [:m1, :n1, :m2, :n2, :h1, :h2, :q]
 - MJO_params: fields are all parameters. Put in dimensional values and output non-dimensional values. (Includes discretization parameters.)
Sets up a bunch of functions (+, -, *, /, maximum, minimum, isnan, isinf) that work on MJO_State and Array{MJO_State,1} variables. 
## time_step.jl
``function f_euler(initial_state:: MJO_State, params::MJO_params, h::Float64, N::Int, every::Int)``
 - This function outputs an Array{MJO_State,1} where it's ith index contains an MJO_State at the every * i + 1 th time-step.
## BOTH mjo_sim_pyplot.jl and mjo_sim_plots_gr.jl 
h: time-step
N: max number of time-steps
every: we only save every "every" time-steps.. 
str: every "every" time-steps of all 6 variables are saved in the folder "../movies/field/str.png"
``function f_euler_contour(
    initial_state:: MJO_State, 
    params::MJO_params, 
    h::Float64, 
    N::Int,
    every::Int,
    str::String
    )``
 - This function saves a png every "every" time-steps until N or until there's a NaN or Inf value in any of the 7 variables. 
 
``function savecontourmaps(evol::Array{MJO_State,1}, str::String; draw::Symbol=:contourf)``
- This function takes an Array{MJO_State,1} (such as the output from f_euler), then saves all of the variables at all times as png's. 

PyPlot version: 
- "draw" is a keyword argument that is by default contourf (filled contour lines). Another possible value is: pcolormesh.
- "draw": imshow will not work, and contour will not work when the 2d array is a constant. 
- PRO: better quality images
- CON: slower.
- CON: can't do much else on the computer while this is running because the figure window keeps popping out!

Plots with gr backend version:
- "draw" is a keyword argument that is by default contour (unfilled, becase filled option looks exactly like heatmap)
- "draw" another possible value is heatmap.
- PRO : faster
- CON : limited quality 
- CON : It crashes pretty often and ends the Julia REPL session completely..
