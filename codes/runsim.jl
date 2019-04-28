include("time_step3.jl")
h_time = 0.0009; #3 min 0.0045 #15 min  
every = 10; #30min #96 * 7; # 1 week
N = every * 8;#every * 52 * 4 #4years
H1 = 0.3;
X=[:T_RC]
x=[gen_params(h_time).T_RC/2]
imex_print(N, every, h_time, "testdom", X=X, x=x, H1=H1, NA=0.4, hov=true) 