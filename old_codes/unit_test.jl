# Tendency Unit Test
# Author: Minah Yang <lucia.yang@colorado.edu>
# Copyright 2018. Not for public distribution.

include("mjo_a.jl")
#= The above file creates params and MJOState structs, the dxdt function, and all other inline functions that are needed for dxdt. It only sets the parameters as params, and sets grid_y and grid_x. 
Note that all tests should check if (all grid points - ghost cells) are working properly. 
=#


using Printf

#######
err = 0.00001
#######



###############################################################################
#=
Test A: Testing Precipitation (C&M) & Radiative Cooling & Divergences == zero. 
If the momenta are zero and heights are constant, and Q is linear from 0 to 1, then the tendency for :
1) m1 is zero,
2) n1 is zero,
3) m2 is zero,
4) n2 is zero,
5) q is just the precipitation,
6) h1 is precipitation and radiative cooling, and
7) h2 is precipitation and radiative cooling.
=#
###############################################################################

println("Test A: Zero momenta, Constant heights, linear Q. This tests for precipiation and radiative cooling while the div_flux terms are turned off.")

state = MJO_State(
          zeros(grid_y, grid_x),                                       # m1
          zeros(grid_y, grid_x),                                       # n1
          zeros(grid_y, grid_x),                                       # m2
          zeros(grid_y, grid_x),                                       # n2
          ones(grid_y, grid_x),                                        # h1
          3.0* ones(grid_y, grid_x),                                   # h2
          repeat(range(10.0^(-16.0), stop=1.1*params.Qs, length=grid_y), 1, grid_x) # q
    );  

outA = deepcopy(state)
dxdt(params, state, outA);
# 1 #
@printf("The max abs value for the error in the tendency for m1 is %.3e. \n", maximum(abs.(outA.m1[2:end-1,:])))

# 2 #
@printf("The max abs value for the error in the tendency for n1 is %.3e. \n", maximum(abs.(outA.n1[2:end-1,:])))

# 3 #
@printf("The max abs value for the error in the tendency for m2 is %.3e. \n", maximum(abs.(outA.m2[2:end-1,:])))

# 4 #
@printf("The max abs value for the error in the tendency for n2 is %.3e. \n", maximum(abs.(outA.n2[2:end,:])))

# 5 #
q_check = (-1.0 .+ params.Qs./(params.DD*state.q)).*(8.0*params.LL)/(params.QQ*86400000.0*params.UU).*(exp.(params.B*state.q/params.Qs).-1.0);
@printf("The max abs discrepancy between P and the tendency for q is %.3e. \n",maximum(abs.( outA.q[2:end-1,:] - q_check[2:end-1,:] )))

# 6 #
rc_check = -params.BQH*(8.0*params.LL)/(params.QQ*86400000.0*params.UU).*(exp.(params.B*state.q/params.Qs).-1.0) + params.Tratio*(state.h2-state.h1)
@printf("The max abs discrepancy between P_RC value and the tendency for h1 is %.3e. \n", maximum(abs.( outA.h1[2:end-1,:] - rc_check[2:end-1,:] )))

# 7 #
@printf("The max abs discrepancy between P_RC value and the tendency for h2 is %.3e. \n", maximum(abs.( outA.h2[2:end-1,:] + rc_check[2:end-1,:] )))




###############################################################################
#= 
Test B: Testing div_flux. 
If (m_1, m_2) are linear in y, (n_1, n_2) are linear in x, Q is (near) 0, (h1,h2) are constant, then the tendency for:
1) m1 should be the Coriolis force + partial_y(v*m_1),
2) n1 should be the Coriolis force + partial_x(u*n_1),
3) m2 should be the Coriolis force + partial_y(v*m_2),
4) n2 should be the Coriolis force + partial_x(u*n_2),
5) q we don't care about?,
6) h1 is 0, and 
7) h2 is 0.
=#
###############################################################################

println("\n\nTest B: Momenta linear in orthogonal directions, Constant heights and (near) Zero Q. This tests the div_flux terms.")


state = MJO_State(
          repeat(range(1,stop=2, length=grid_y), 1, grid_x), # m1
          repeat(range(3,stop=4, length=grid_x)', grid_y,1), # n1
          repeat(range(1,stop=2, length=grid_y), 1, grid_x), # m2
          repeat(range(3,stop=4, length=grid_x)',grid_y, 1), # n2
          ones(grid_y, grid_x),                     # h1
          ones(grid_y, grid_x),                     # h2
          zeros(grid_y,grid_x)                      # q
    );                                          
outB = deepcopy(state)
dxdt(params, state, outB);

function check_div_flux(m, n, h, q, delt_x, delt_y)
  m = [m[:,end] m m[:,1]];
  n = [n[:,end] n n[:,1]];
  h = [h[:,end] h h[:,1]];
  q = [q[:,end] q q[:,1]];

  x = .25*delt_x*(
            +(((m[2:end-1,2:end-1]./h[2:end-1,2:end-1]) + (m[2:end-1,3:end]./h[2:end-1,3:end])).*(q[2:end-1,2:end-1]+q[2:end-1,3:end]))
            - (((m[2:end-1,1:end-2]./h[2:end-1,1:end-2]) + (m[2:end-1,2:end-1]./h[2:end-1,2:end-1])).*(q[2:end-1,1:end-2]+q[2:end-1,2:end-1])))
  y = .25*delt_y*(
  (((n[2:end-1,2:end-1]./h[2:end-1,2:end-1]) + (n[3:end,2:end-1]./h[3:end,2:end-1])).*(q[2:end-1,2:end-1]+q[3:end,2:end-1]))
  -(((n[1:end-2,2:end-1]./h[1:end-2,2:end-1]) + (n[2:end-1,2:end-1]./h[2:end-1,2:end-1])).*(q[1:end-2,2:end-1]+q[2:end-1,2:end-1])))
  return x+y
end

# 1 #
q = deepcopy(state.m1)
m1_check = (params.Ro*state.n1[2:end-1,:].*params.y[2:end-1]) - check_div_flux(state.m1, state.n1, state.h1, q, params.delt_x, params.delt_y)
@printf("The max abs discrepancy between -n1/Ro and tendency for m1 is %.3e.\n", maximum(abs.(m1_check -outB.m1[2:end-1,:])))

# 2 #
q = deepcopy(state.n1)
n1_check = (-params.Ro*state.m1[2:end-1,:].*params.y[2:end-1]) - check_div_flux(state.m1, state.n1, state.h1, q, params.delt_x, params.delt_y)
@printf("The max abs discrepancy between +m1/Ro and tendency for n1 is %.3e.\n", maximum(abs.(n1_check -outB.n1[2:end-1,:])))

# 3 #
q = deepcopy(state.m2)
m2_check = (params.Ro*state.n2[2:end-1,:].*params.y[2:end-1]) - check_div_flux(state.m2, state.n2, state.h2, q, params.delt_x, params.delt_y)
@printf("The max abs discrepancy between -n2/Ro and tendency for m2 is %.3e.\n", maximum(abs.(m2_check -outB.m2[2:end-1,:])))

# 4 #
q = deepcopy(state.n2)
n2_check = (-params.Ro*state.m2[2:end-1,:].*params.y[2:end-1]) - check_div_flux(state.m2, state.n2, state.h2, q, params.delt_x, params.delt_y)
@printf("The max abs discrepancy between +m2/Ro and tendency for n2 is %.3e.\n", maximum(abs.(n2_check -outB.n2[2:end-1,:])))

# 5 #
println("Since we've set q=0, we can't test it.")

# 6 #
@printf("The max abs value for the error in the tendency for h1 is %.3e. \n", maximum(abs.(outB.h1[2:end-1,:])))

# 7 #
@printf("The max abs value for the error in the tendency for h2 is %.3e. \n", maximum(abs.(outB.h2[2:end-1,:])))


###############################################################################
#= 
Test C: Still testing div_flux. 
Let the heights be constant and Q = 0. If m1, m2 are sin(2pix/LX) and the other momenta are zero, then the tendency for:
1) m1 should be proprotional to 2sin(2pix/LX)cos(2pix/LX),
2) n1 should be the Coriolis Force,
3) m2 should be proprotional to 2sin(2pix/LX)cos(2pix/LX),
4) n2 is should be the Coriolis Force,
5) q we don't care about?,
6) h1 should be proportional to cos(2pix/LX), and 
7) h2 should be proportional to cos(2pix/LX).
=#
###############################################################################


println("\n\nTest C: Zonal momenta sinusoidal along latitude, Constant heights and (near) Zero Q. This tests the div_flux terms acting as a simple derivative.")

x = range(0, stop=params.lon[end]*pi/180, length=grid_x)
# m1 will be sin(x/T) where the period T is given by the nondimensionalized length of the zonal domain.
# T = 2*pi* RE/ LL
# i.e. T = RE/ LL, since x ranges from 0 to 2*pi-delt_x. 
state = MJO_State(
          repeat(sin.(x)', grid_y,1), # m1
          zeros(grid_y, grid_x),      # n1
          repeat(sin.(x)', grid_y,1), # m2
          zeros(grid_y, grid_x),      # n2
          ones(grid_y, grid_x),       # h1
          ones(grid_y, grid_x),       # h2
          zeros(grid_y, grid_x)       # q
    );                                          
outC = deepcopy(state)
dxdt(params, state, outC);
TC = params.RE/params.LL;

# 1 #
ind = []
for i = 1 : grid_x
    if 1-abs(-TC*outC.m1[2,i]/(sin(x[i])*cos(x[i]))) > err
        push!(ind,i)
    end
end
if length(ind)>0
  println("m1: The inaccurate indices are:",ind, "and their respective locations in [0,2pi] are:", ind/grid_x)
else
  @printf("m1: All entries were within %.3e percent relative error.\n", 100.0*err)
end

# 2 #
check_n1 = -params.Ro*state.m1[2:end-1,:].*params.y[2:end-1]
@printf("The max abs value for the error in the tendency for n1 is %.3e. \n", maximum(abs.(check_n1 - outC.n1[2:end-1,:])))

# 3 #
ind = []
for i = 1 : grid_x
    if 1-abs(-TC*outC.m2[2,i]/(sin(x[i])*cos(x[i]))) > err
        push!(ind,i)
    end
end
if length(ind)>0
  println("m2: The inaccurate indices are:",ind, "and their respective locations in [0,2pi] are:", ind/grid_x)
else
  @printf("m2: All entries were within %.3e percent relative error.\n", 100.0*err)
end

# 4 #
check_n2 = -params.Ro*state.m2[2:end-1,:].*params.y[2:end-1]
@printf("The max abs value for the error in the tendency for n2 is %.3e. \n", maximum(abs.(check_n2 - outC.n2[2:end-1,:])))

# 5 #
println("Since we've set q=0, we can't test it.")

# 6 #
ind = []
for i = 1 : grid_x
    if 1-abs(-TC*outC.h1[2,i]/cos(x[i])) > err
        push!(ind,i)
    end
end
if length(ind)>0
  println("h1: The inaccurate indices are:",ind, "and their respective locations in [0,2pi] are:", ind/grid_x)
else
  @printf("h1: All entries were within %.3e percent relative error.\n", 100.0*err)
end

# 7 #
ind = []
for i = 1 : grid_x
    if 1-abs(-TC*outC.h2[2,i]/cos(x[i])) > err
        push!(ind,i)
    end
end
if length(ind)>0
  println("h2: The inaccurate indices are:",ind, "and their respective locations in [0,2pi] are:", ind/grid_x)
else
  @printf("h2: All entries were within %.3e percent relative error.\n", 100.0*err)
end

###############################################################################
#= 
Test D: Still testing div_flux. 
Let the heights be constant and Q = 0. If n1, n2 are sin(2piy/LY) and the other momenta are zero, then the tendency for:
1) m1 is just the Coriolis force,
2) n1 should be proprotional to sin(2pix/LX)cos(2pix/LX),
3) m2 is just the Coriolis force,
4) n2 should be proprotional to sin(2pix/LX)cos(2pix/LX),
5) q we don't care about?,
6) h1 should be proportional to cos(2pix/LX), and 
7) h2 should be proportional to cos(2pix/LX).
=#
###############################################################################


println("\n\nTest D: Meridional momenta sinusoidal longitude, Constant heights and (near) Zero Q. This tests the div_flux terms acting as a simple derivative.")

y = range(0, stop = 2*pi, length = grid_y);
# n1 will be sin(y/T) where the period T is given by the nondimensionalized length of the meridional domain. 
# T = (20.125-(-20.125))*pi/180 * RE/LL
# i.e. T = (20.125-(-20.125))*pi/180 * RE/LL, since y ranges from 0 to 2*pi. 
TD = (20.125-(-20.125))*pi/180 * params.RE/params.LL;
state = MJO_State(
          zeros(grid_y, grid_x),      # m1
          repeat(sin.(y), 1, grid_x), # n1
          zeros(grid_y, grid_x),      # m2
          repeat(sin.(y), 1, grid_x), # n2
          ones(grid_y, grid_x),       # h1
          ones(grid_y, grid_x),       # h2
          zeros(grid_y, grid_x)       # q
    );                                          
outD = deepcopy(state);
dxdt(params, state, outD);

# 1 #
check_m1 = params.Ro*state.n1[2:end-1,:].*params.y[2:end-1]
@printf("The max abs value for the error in the tendency for m1 is %.3e. \n", maximum(abs.(check_m1 - outD.m1[2:end-1,:])))

# 2 #
ind = []
for i = 2 : grid_y-1
    if 1-abs(TD*outD.n1[i,2]/(sin(y[i])*cos(y[i]))) > err
        push!(ind,i)
    end
end
if length(ind)>0
  println("n1: The inaccurate indices are:",ind, "and their respective locations in [0,2pi] are:", ind/grid_x)
else
  @printf("n1: All entries were within %.3e percent relative error.\n", 100.0*err)
end


# 3 #
check_m2 = params.Ro*state.n2[2:end-1,:].*params.y[2:end-1]
@printf("The max abs value for the error in the tendency for m2 is %.3e. \n", maximum(abs.(check_m2 - outD.m2[2:end-1,:])))

# 4 #
ind = []
for i = 2 : grid_y-1
    if 1-abs(TD*outD.n2[i,2]/(sin(y[i])*cos(y[i]))) > err
        push!(ind,i)
    end
end
if length(ind)>0
  println("n2: The inaccurate indices are:",ind, "and their respective locations in [0,2pi] are:", ind/grid_x)
else
  @printf("n2: All entries were within %.3e percent relative error.\n", 100.0*err)
end

# 5 #
println("Since we've set q=0, we can't test it.")

# 6 #
ind = []
for i = 1 : grid_x
    if 1-abs(-TD*outD.h1[2,i]/cos(x[i])) > err
        push!(ind,i)
    end
end
if length(ind)>0
  println("h1: The inaccurate indices are:",ind, "and their respective locations in [0,2pi] are:", ind/grid_x)
else
  @printf("h1: All entries were within %.3e percent relative error.\n", 100.0*err)
end

# 7 #
ind = []
for i = 1 : grid_x
    if 1-abs(-TD*outD.h2[2,i]/cos(x[i])) > err
        push!(ind,i)
    end
end
if length(ind)>0
  println("h2: The inaccurate indices are:",ind, "and their respective locations in [0,2pi] are:", ind/grid_x)
else
  @printf("h2: All entries were within %.3e percent relative error.\n", 100.0*err)
end



###############################################################################
#= 
Test E: Testing Gradient. 
Let the heights sinusoidal in x, and everything else 0. Then the tendency for:
1) m1 should be linear,
2) n1 should be 0,
3) m2 should be linear,
4) n2 should be 0,
5) q we don't care about?,
6) h1 should be 0,
7) h2 should be 0.
=#
###############################################################################

println("\n\nTest E: Momenta zero, heights sinusoidal along x, and Zero Q. This tests the gradient terms acting as a simple derivative.")

x = range(0, stop=params.lon[end]*pi/180, length=grid_x)
# m1 will be sin(x/T) where the period T is given by the nondimensionalized length of the zonal domain.
# T = 2*pi* RE/ LL
# i.e. T = RE/ LL, since x ranges from 0 to 2*pi-delt_x. 
TE = params.LL/params.RE
state = MJO_State(
          zeros(grid_y, grid_x),                    # m1
          zeros(grid_y, grid_x),                    # n1
          zeros(grid_y, grid_x),                    # m2
          zeros(grid_y, grid_x),                    # n2
          2 .+repeat(sin.(x)', grid_y,1), # h1
          2 .+repeat(sin.(x)', grid_y,1), # h2
          zeros(grid_y, grid_x),                    # q
    );                                          
outE = deepcopy(state)
dxdt(params, state, outE);


# 1 #
check_m1 = -params.Fr*2.0*state.h1.*repeat(cos.(x)', grid_y,1)*TE
@printf("The max abs value for the error in the tendency for m1 is %.3e. \n", maximum(abs.(check_m1[2:end-1,:] - outE.m1[2:end-1,:])))

# 2 #
@printf("The max abs value for the error in the tendency for n1 is %.3e. \n", maximum(abs.(outE.n1[2:end-1,:])))

# 3 #
check_m2 = -params.Fr*(1.0+params.AA)*state.h1.*repeat(cos.(x)', grid_y,1)*TE
@printf("The max abs value for the error in the tendency for m2 is %.3e. \n", maximum(abs.(check_m2[2:end-1,:] - outE.m2[2:end-1,:])))

# 4 #
@printf("The max abs value for the error in the tendency for n2 is %.3e. \n", maximum(abs.(outE.n2[2:end-1,:])))

# 5 #
println("Since we've set q=0, we can't test it.")

# 6 #
@printf("The max abs value for the error in the tendency for h1 is %.3e. \n", maximum(abs.(outE.h1[2:end-1,:])))

# 7 #
@printf("The max abs value for the error in the tendency for h2 is %.3e. \n", maximum(abs.(outE.h2[2:end-1,:])))

###############################################################################
#= 
Test F: Testing Gradient. 
Let the heights linear in x, and everything else 0. Then the tendency for:
1) m1 should be linear,
2) n1 should be 0,
3) m2 should be linear,
4) n2 should be 0,
5) q we don't care about?,
6) h1 should be 0,
7) h2 should be 0.
=#
###############################################################################

println("\n\nTest F: Momenta zero, heights linear along x, and Zero Q. This tests the gradient terms acting as a simple derivative.")

# n1 will be sin(y/T) where the period T is given by the nondimensionalized length of the meridional domain. 
# T = (20.125-(-20.125))*pi/180 * RE/LL
# i.e. T = (20.125-(-20.125))*pi/180 * RE/LL, since y ranges from 0 to 2*pi. 
TF = 1.0/( (20.125-(-20.125))*pi/180 * params.RE/params.LL );

state = MJO_State(
          zeros(grid_y, grid_x),                    # m1
          zeros(grid_y, grid_x),                    # n1
          zeros(grid_y, grid_x),                    # m2
          zeros(grid_y, grid_x),                    # n2
          repeat(range(1,stop=2,length=grid_y), 1, grid_x), # h1
          repeat(range(1,stop=2,length=grid_y), 1, grid_x), # h2
          zeros(grid_y, grid_x),                    # q
    );                                          
outF = deepcopy(state)
dxdt(params, state, outF);

# 1 #
@printf("The max abs value for the error in the tendency for m1 is %.3e. \n", maximum(abs.(outF.m1[2:end-1,:])))

# 2 #
check_n1 = -params.Fr*2.0*state.h1*TF
@printf("The max abs value for the error in the tendency for n1 is %.3e. \n", maximum(abs.(check_n1[2:end-1,:] - outF.n1[2:end-1,:])))

# 3 #
@printf("The max abs value for the error in the tendency for m2 is %.3e. \n", maximum(abs.(outF.m2[2:end-1,:])))

# 4 #
check_n2 = -params.Fr*(1.0+params.AA)*state.h1*TF
@printf("The max abs value for the error in the tendency for n2 is %.3e. \n", maximum(abs.(check_n2[2:end-1,:] - outF.n2[2:end-1,:])))

# 5 #
println("Since we've set q=0, we can't test it.")

# 6 #
@printf("The max abs value for the error in the tendency for h1 is %.3e. \n", maximum(abs.(outF.h1[2:end-1,:])))

# 7 #
@printf("The max abs value for the error in the tendency for h2 is %.3e. \n", maximum(abs.(outF.h2[2:end-1,:])))