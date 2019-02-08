include("time_step3.jl")
h_time = 0.0045; #15 iin
every = 96;      #24 hrs
N = every*365*2;  #2 years
imex_print(N, every, h_time, "15min-24hr", step=3)
