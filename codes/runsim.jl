include("time_step3.jl")
H1=0.3; 
T_RC = 8*24*3600.0; #8 days
h_time = 0.0009; #3 min
PP=12500.0;
lat_range=[-30.0, 30.0];
deg = 0.25;
DD = 1.5;
bb = 0.004*h_time; #(KK*h_time)
params=gen_params(H1=H1, DD=DD, T_RC=T_RC, PP=PP, lat_range=lat_range, deg=deg, bb=bb, h_time=h_time)
#added moisture source in q tendency

every = 5; #1hour
N = every * 8; #* 24 * 365; #1year
NA = 0.4; 
fr = 0.5;
msource=[0.64, 3.0]; #amplitude, # of stdev by 23.4 deg. 
hovtxtfile = "asdf" #Name of hovmoller info txt file
hovpngname = "asdf"  #Name of hovmoller png file

IC = genInitSr(name="mresend")
#=evol = imex(N, every, h_time,
	bb=bb, params=params, NA=NA, fr=fr, msource=msource);
x=0;=#
imex_print(N, every, h_time, hovtxtfile, 
	bb=bb, params=params, IC=IC, NA=NA, fr=fr, msource=msource,
	hov=true, everyH=every)
hovmollertxt(hovtxtfile, hovpngname,T=60)
