include("time_step3.jl")
H1=0.3; 
T_RC = 8*24*3600.0; #8 days
h_time = 0.0009; #3 min
PP=22500.0;
params=gen_params(H1=H1, T_RC=T_RC, PP=PP, h_time=h_time)

every = 20; #30min #96 * 7; # 1 week
N = every * 8;#every * 52 * 4 #4years


bb = 0.005;
NA = 0.4; 
hovtxtfile = "testfric" #Name of hovmoller info txt file
hovpngname = "hovfric"  #Name of hovmoller png file

imex_print(N, every, h_time, "testfric", bb=bb, params=params, NA=NA, hov=true, everyH=every) 
hovmollertxt("testfric", "hovfric", T=120)