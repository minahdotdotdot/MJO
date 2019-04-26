include("time_step3.jl")
h_time = 0.0045 #15 min 0.0009; #3 min  
every = 96 * 7; # 1 week
N = every * 52 * 4 #4years
H1 = 0.03;
NA=range(0.03, stop=0.2, length=9)
imex_print(N, every, h_time, "15min-12hrH",H1=H1, NA=NA[1]) 
