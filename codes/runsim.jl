include("time_step3.jl")
h_time = 0.0045 #15 min 0.0009; #3 min  
every = 48 #12 hours 2.4 hrs
N = every*365*2 #a day.#every#*365*2;  #2 years
X = [:B, :PP];                       #parameters to tweak
H1 = range(1.8, stop=1.2, length=9);
imex_print(N, every, h_time, "15min-12hrH",H1=1.8)                                                   
