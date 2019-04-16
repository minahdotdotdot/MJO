include("time_step3.jl")
h_time = 0.0045 #15 min 0.0009; #3 min  
every = 48 #12 hours 2.4 hrs
N = every*365*2 #a day.#every#*365*2;  #2 years
X = [:B, :PP];                       #parameters to tweak
B  = range(11.0,    stop=11.9,    length=3);
PP = range(15000.0, stop=22000.0, length=3);
x = [B[1], PP[1]];
imex_print(N, every, h_time, "15min-12hrH", step=3, X=X, x=x, H1=1.8)                                                   
