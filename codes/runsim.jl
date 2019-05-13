include("time_step3.jl")
H1=0.3; 
T_RC = 8*24*3600.0; #8 days
h_time = 0.0009; #3 min
PP=12500.0;
lat_range=[-30.0, 30.0];
deg = 0.0625;
DD = 1.5;
params=gen_params(H1=H1, DD=DD, T_RC=T_RC, PP=PP, lat_range=lat_range, deg=deg, h_time=h_time)
#added moisture source in q tendency
every = 5; #1hour
N = every * 8; #* 24 * 365; #1year


bb = 0.005;
NA = 0.4; 
fr = 0.5;
hovtxtfile = "test" #Name of hovmoller info txt file
hovpngname = "testhov"  #Name of hovmoller png file

imex_print(N, every, h_time, hovtxtfile, bb=bb, params=params, NA=NA, fr=fr, hov=true, everyH=every) 
hovmollertxt(hovtxtfile, hovpngname,T=60)
