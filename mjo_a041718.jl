struct MJO_params
	LL    :: Float64 #\hat{L}     ~ 10^6              (m) 
	HH    :: Float64 #H           ~ 5000              (m) 
	UU    :: Float64 #U           ~ 5                 (m/s)
	QQ    :: Float64 #\hat{Q}     ~ 0.05              (m)
	T_RC  :: Float64 #T_{RC}      ~ 16 days = 1382400 (s)
	T_Q   :: Float64 #\tau_q      ~	????????          (s)

	g     :: Float64 # gravity on earth : 9.80665     (m/s^2)
	RE    :: Float64 # radius of earth  : 6,371,000   (m)

	AA    :: Float64 # \alpha \approx 1.0183548918336027
	BB    :: Float64 # \beta \approx 750.
	DD    :: Float64 # CM \delta \approx 1.1
	Qs    :: Float64 # moisture saturation 56.2 to 70. is reasonable
	B     :: Float64 # CM "b" \approx 11.4

	Ro    :: Float64 # \frac{1}{Ro*y} = 4*pi*L^2/(U*day in seconds * RE) 
	Fr    :: Float64 # 1/Fr^2 = gH/U^2
	BQH   :: Float64 #\frac{\hat{Q}}{H}
	Tratio:: Float64 # :L / (U*T_RC)

	MJO_params(LL, HH, UU, QQ, T_RC, T_Q, AA, BB, DD, Qs, B) = {
	new(copy(LL), copy(HH), copy(UU), copy(QQ), copy(T_RC), copy(T_Q),
		copy(g), copy(RE), 
		copy(AA), copy(BB), copy(DD), copy(Qs), copy(B),
		4.*pi*LL^2/(3600.*24.*UU*RE), g*HH/UU^2, BB*QQ/HH, L/(U*T_RC))
	}
end

struct MJO_State
	delt_x :: Float64
	delt_y :: Float64
	lon    :: Array{Float64, 1}
	lat    :: Array{Float64, 1}

	y 	   :: Array{Float64, 1}

 	m1     :: Array{Float64, 2}
	n1     :: Array{Float64, 2}
	m2     :: Array{Float64, 2}
	n2     :: Array{Float64, 2}
	h1     :: Array{Float64, 2}
	h2     :: Array{Float64, 2}
	q      :: Array{Float64, 2}

	MJO_state(delt_x, delt_y, lon, lat, m1, n1, m2, n2, h1, h2, q) = {
	new(copy(delt_x), copy(delt_y), copy(lat), copy(lon), copy(y),
		copy(m1), copy(n1), copy(m2), copy(n2), copy(h1), copy(h2), copy(q))
	}
end

params = MJO_params(10.^6, 5000., 5., 0.05, 1382400., 1., 1.0184, 750., 1.1, 58., 11.4) #FIND tau_q (T_Q)
state = {MJO_state(.25*pi/180*params.RE/params.LL, .25*pi/180*params.RE/params.LL, linspace(0., 360.-.25, 1440), linspace(-20.+.125, 20.-.125, 160),
	pi/180*params.RE * linspace(0., 360.-.25, 1440)/params.LL, 
	zeros(1440, 160), zeros(1440, 160), zeros(1440, 160), zeros(1440, 160), zeros(1440, 160), zeros(1440, 160), zeros(1440, 160))
}


##Craig and Mack Precipitation 
@inline function P{T}(LL::T, UU::T, QQ::T, B::T, Qs:: T, q::T, T_Q::T)# add boolean/string to choose between {C&M ,Betts-Miller}? 
	value = (8./QQ)/(86400000.*UU/LL)*(exp(B*q/Qs)-1.)
	#=
	value = 0
	if q > Qs
		value = ((LL/UU)/T_Q)*(q-Qs)
	end
	=#
	return value

## Precipitation and Radiative Cooling together (RC never appears on its own.)
@inline function P_RC{T} (LL::T, UU::T, QQ::T, B::T, Qs:: T, q::T, T_Q::T, BQH::T, Tratio::T, hh2::T, hh1::T )
	value = P(LL,UU,QQ,B,Qs,q,T_Q)
	if hh2 - hh1 > T(0)
		value = value - Tratio * (hh2-hh1)
	return value
end

@inline function div_flux{T} (m:: Array{Float64, 2}, n:: Array{Float64, 2}, h:: Array{Float64, 2}, q:: Array{Float64, 2}, ii::Int64, jj::Int64, delt_x::Float64, delt_y::Float64)
	#Account for zonal periodicity.
	if ii == 1
		iii = 1441 # to satisfy iii-1 = 1440
	end
	if ii == 1440
		iiii = 0 #to satisfy iiii+1 = 1
	end
	value = 0 

	if jj == 1
		### LOWER BOUNDARY ###
		value = {
		.25/delt_x*(
			+(m[jj,ii]/h[jj,ii] + m[jj,ii+1]/h[jj,ii+1])*(q[jj,ii]+q[jj,ii+1]) 
			- ((m[jj,ii-1]/h[jj,ii-1] + m[jj,ii]/h[jj,ii]))(q[jj,ii-1]+q[jj,ii])
		)+.25/delt_y*(
			+(n[jj,ii]/h[jj,ii] + n[jj+1,ii]/h[jj+1,ii])*(q[jj,ii]+q[jj+1,ii])
		)
			}
	else if jj == 360
		### UPPER BOUNDARY ###
		value = {
		.25/delt_x*(
			+(m[jj,ii]/h[jj,ii] + m[jj,ii+1]/h[jj,ii+1])*(q[jj,ii]+q[jj,ii+1]) 
			- ((m[jj,ii-1]/h[jj,ii-1] + m[jj,ii]/h[jj,ii]))(q[jj,ii-1]+q[jj,ii])
		)+.25/delt_y*(
			-(n[jj-1,ii]/h[jj-1,ii] + n[jj,ii]/h[jj,ii])*(q[jj-1,ii]+q[jj,ii])
		)
			}
	else 
		### INTERIOR ###
		value = {
		.25/delt_x*(
			+(m[jj,ii]/h[jj,ii] + m[jj,ii+1]/h[jj,ii+1])*(q[jj,ii]+q[jj,ii+1]) 
			- ((m[jj,ii-1]/h[jj,ii-1] + m[jj,ii]/h[jj,ii]))(q[jj,ii-1]+q[jj,ii])
		)+.25/delt_y*(
			+(n[jj,ii]/h[jj,ii] + n[jj+1,ii]/h[jj+1,ii])*(q[jj,ii]+q[jj+1,ii])
			-(n[jj-1,ii]/h[jj-1,ii] + n[jj,ii]/h[jj,ii])*(q[jj-1,ii]+q[jj,ii])
		)
			}
	return value

function dxdt(params:: MJO_params, state::MJO_State, out::MJO_State)
	delt_x = state.delt_x  ::Float64
	delt_y = state.delt_y  ::Float64
	LL     = params.LL     ::Float64
	UU     = params.UU     ::Float64
	QQ     = params.QQ     ::Float64
	B      = params.B      ::Float64
	Qs     = params.T      ::Float64
	BQH    = params.BQH    ::Float64
	Tratio = params.Tratio ::Float64
	AA     = params.AA     ::Float64
	DD     = params.DD     ::Float64
	T_Q    = params.T_Q    ::Float64

	for ii = 1 : length(MJO_state.lon)	
		iii = ii  #Use for ii-1
		iiii = ii #Use for ii+1
		#Periodic boundary conditions in lon: lon[1] = 0 = 360, lon[end] = long[1440] = 359.75
		if ii == 1
			iii = 1441 # to satisfy iii-1 = 1440
		end
		if ii == 1440
			iiii = 0 #to satisfy iiii+1 = 1
		end

		value_P_RC = P_RC(LL, UU, QQ, B, Qs, state.q[jj,ii], T_Q, BQH, Tratio, state.h2[jj,ii], state.h1[jj,ii])

		
		for jj = 1: length(MJO_state.y) 
			### MOMENTUM
			out.m1[jj,ii] = {
			+ div_flux(state.m1, state.n1, state.h1, state.m1, ii, jj, delt_x, delt_y)
			+ params.Ro*state.y[jj,ii]*state.n1[jj,ii] #=1/Ro*n1=#
			- params.Fr*( state.h1[jj,ii] * .5*delt_x*(
				+(state.h1[jj,ii+1] + state.h2[jj,ii+1])
				-(state.h1[jj,ii-1] + state.h2[jj,ii-1])
				)#=h_1∂x(h_1+h_2)=#
			) 
			- state.m1[jj,ii]/state.h1[jj,ii]*value_P_RC 
			}

			out.n1[jj,ii] = {
			+ div_flux(state.m1, state.n1, state.h1, state.n1, ii, jj, delt_x, delt_y)
			- params.Ro*state.y[jj,ii]*state.m1[jj,ii] #=-1/Ro*m1=#
			- params.Fr*( state.h1[jj,ii] * .5*delt_y*(
				+(state.h1[jj+1,ii] + state.h2[jj+1,ii])
				-(state.h1[jj-1,ii] + state.h2[jj-1,ii])
				) #=h_1∂y(h_1+h_2)=#
			) 
			- state.n1[jj,ii]/state.h1[jj,ii]*value_P_RC 
			}

			out.m2[jj,ii] = {
			+ div_flux(state.m2, state.n2, state.h2, state.m2, ii, jj, delt_x, delt_y)
			+ params.Ro*state.y[jj,ii]*state.n2[jj,ii] #=1/Ro*n2=#
			- params.Fr*( state.h1[jj,ii] * .5*delt_x*(
				+(state.h1[jj,ii+1] + AA*state.h2[jj,ii+1])
				-(state.h1[jj,ii-1] + AA*state.h2[jj,ii-1])
				) #=h_1∂x(h_1+\alpha* h_2)=#
			)
			+ state.m1[jj,ii]/state.h1[jj,ii]*value_P_RC 
			}

			out.n2[jj,ii] = {
			+ div_flux(state.m2, state.n2, state.h2, state.n2, ii, jj, delt_x, delt_y)
			- params.Ro*state.y[jj,ii]*state.m2[jj,ii] #=-1/Ro*m2=#
			- params.Fr*( state.h1[jj,ii] * .5*delt_y*(
				+(state.h1[jj+1,ii] + AA*state.h2[jj+1,ii])
				-(state.h1[jj-1,ii] + AA*state.h2[jj-1,ii])
				) #=h_1∂y(h_1+\alpha*h_2)=#
			)
			+ state.n1[jj,ii]/state.h1[jj,ii]*value_P_RC 
			}

			### MASS

			if jj == 1

				######################	
				### LOWER BOUNDARY ###   jj-1 shouldn't exist in this section.
				######################

				out.h1[jj,ii] = {
				-.5/delt_x * (state.m1[jj,iiii+1] - state.m1[jj,iii-1]) #=∂x m_1=#
				-.5/delt_y * ((state.n1[jj+1,ii] + state.n1[jj,ii])-(0.)) #=∂y n_1=#
				- value_P_RC
				}

				out.h2[jj,ii] = {
				-.5/delt_x * (state.m2[jj,iiii+1] - state.m2[jj,iii-1]) #=∂x m_2=#
				-.5/delt_y * ((state.n2[jj+1,ii] + state.n2[jj,ii])-(0.)) #=∂y n_2=#
				+ value_P_RC
				}
			else if jj == length(state.y)	

				######################	
				### UPPER BOUNDARY ###
				######################

				out.h1[jj,ii] = {
				-.5/delt_x * (state.m1[jj,iiii+1] - state.m1[jj,iii-1]) #=∂x m_1=#
				-.5/delt_y * ((0.)-(state.n1[jj-1,ii] + state.n1[jj,ii])) #=∂y n_1=#
				- value_P_RC
				}

				out.h2[jj,ii] = {
				-.5/delt_x * (state.m2[jj,iiii+1] - state.m2[jj,iii-1]) #=∂x m_2=#
				-.5/delt_y * ((0.)-(state.n2[jj-1,ii] + state.n2[jj,ii]) #=∂y n_2=#
				+ value_P_RC
				}
			else 

				######################	
				###### INTERIOR ######
				######################

				out.h1[jj,ii] = {
				-.5/delt_x * (state.m1[jj,iiii+1] - state.m1[jj,iii-1]) #=∂x m_1=#
				-.5/delt_y * (state.n1[jj+1,ii] - state.n1[jj-1,ii]) #=∂y n_1=#
				- value_P_RC
				}

				out.h2[jj,ii] = {
				-.5/delt_x * (state.m2[jj,iiii+1] - state.m2[jj,iii-1]) #=∂x m_2=#
				-.5/delt_y * (state.n2[jj+1,ii] - state.n2[jj-1,ii]) #=∂y n_2=#
				+ value_P_RC
				}
			end

			### MOISTURE
			out.q[jj,ii]  = {
			+ div_flux(state.m1, state.n1, state.h1, state.q, ii, jj, delt_x, delt_y)
			+(-1+1/(DD*state.q[jj,ii]/Qs))*P(LL,UU,QQ,B,Qs,state.q[jj,ii], T_Q) #=\hat{P}(Q)=#
			}
	end
	return out
end
