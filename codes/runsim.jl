include("time_step3.jl")
h_time = 0.0045; #15 iin
every = 96;      #24 hrs
N = every*365*2;  #2 years
X = [:B, :PP]; 
x = [11.3, 19500.0]
imex_print(N, every, h_time, "15min-24hr", step=3, X, x)
