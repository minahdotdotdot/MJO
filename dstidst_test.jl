using FFTW, LinearAlgebra
N = 14400

# Normal DST
x = rand(N); 
xhat = FFTW.r2r(x, FFTW.RODFT10); 
xx = .5/N*FFTW.r2r(xhat, FFTW.RODFT01); 
print("normal: ", norm(x-xx,2), "\n")

# Dropping last coeff.
xhat = vcat(0, FFTW.r2r(x, FFTW.RODFT10)[1:end-1]); 
xx = .5/N*FFTW.r2r(vcat(
	xhat[2:end],
	0
	), FFTW.RODFT01); 
print("dropped last coeff: ", norm(x-xx,2), "\n")

# Below block of code just shows difference in normal idst(dst) vs. truncated idst(dst)
#=
Trials = 500
err = zeros(2,Trials)
for i = 1 : Trials
	x = rand(N); 
	xhat = FFTW.r2r(x, FFTW.RODFT10); 
	xx = .5/N*FFTW.r2r(xhat, FFTW.RODFT01); 
	err[1,i] = norm(x-xx,2)

	xhat = vcat(0, FFTW.r2r(x, FFTW.RODFT10)[1:end-1]); 
	xx = .5/N*FFTW.r2r(vcat(
		xhat[2:end],
		0
		), FFTW.RODFT01);
	err[2,i] = norm(x-xx,2)
end
using PyPlot
fig, ax = subplots(); 
ax[:set_yscale]("log")
scatter(range(1, stop= Trials, length=Trials), err[1,:], label="normal")
scatter(range(1, stop= Trials, length=Trials), err[2,:], label="trucated")
legend()
=#

#=function dcsft(state::MJO_State, statehat::MJO_State_im,
    fields=fieldnames(MJO_State)[1:6], grid_x::Int=1440)
    for qq in fields
        if qq in fields[[1,3,5,6]]   # m1, m2, h1, h2 
            # DCT then RFFT
            getproperty(statehat,qq)[:,:] = FFTW.rfft(
                vcat(
                FFTW.r2r(getproperty(state,qq), FFTW.REDFT10, 1),
                zeros(1,grid_x)
                ), 2)
        else                        # n1, n2
            # DST then RFFT
            getproperty(statehat,qq)[:,:] = FFTW.rfft(
                vcat(
                    zeros(1,grid_x),
                    FFTW.r2r(getproperty(state,qq), FFTW.RODFT10, 1)
                    ), 2)
        end
    end
    return statehat
end=#

#=
function idcsft(state::MJO_State, statehat::MJO_State_im,
    fields=fieldnames(MJO_State)[1:6], grid_x::Int=1440)
    for qq in fields
        if qq in fields[[1,3,5,6]]   # m1, m2, h1, h2 
            # iRFFT then iDCT
            getproperty(state,qq)[:,:] = .5/grid_y*FFTW.r2r(
                real(
                    FFTW.irfft(getproperty(statehat,qq), grid_x, 2)
                    )[1:end-1,:],
                FFTW.REDFT01, 1)
        else                        # n1, n2
            # iRFFT then iDST
            getproperty(state,qq)[:,:] = .5/grid_y*FFTW.r2r(
                real(
                    FFTW.irfft(getproperty(statehat,qq), grid_x, 2)
                    )[2:end,:],
                FFTW.RODFT01, 1)
        end
    end
    return state
end=#

function dcsft(state::MJO_State, statehat::MJO_State_im,
    fields=fieldnames(MJO_State)[1:6], grid_x::Int=1440)
    for qq in fields
        if qq in fields[[1,3,5,6]]   # m1, m2, h1, h2 
            # DCT then RFFT
            getproperty(statehat,qq)[:,:] = FFTW.rfft(
                vcat(
                    FFTW.r2r(getproperty(state,qq)[2:end-1,:], FFTW.REDFT10, 1),
                    zeros(1,grid_x)
                    ), 2)
        else                        # n1, n2
            # DST then RFFT
            getproperty(statehat,qq)[:,:] = FFTW.rfft(
                vcat(
                    zeros(1,grid_x),
                    FFTW.r2r(getproperty(state,qq)[2:end-1,:], FFTW.RODFT10, 1)#[1:end-1, :]
                    ), 2)
        end
    end
    return statehat
end

function idcsft(state::MJO_State, statehat::MJO_State_im,
    fields=fieldnames(MJO_State)[1:6], grid_x::Int=1440)
    for qq in fields
        if qq in fields[[1,3,5,6]]   # m1, m2, h1, h2 
            # iRFFT then iDCT
            getproperty(state,qq)[2:end-1,:] = .5/(grid_y-2)*FFTW.r2r(
                real(
                    FFTW.irfft(getproperty(statehat,qq), grid_x, 2)
                    )[1:end-1,:],
                FFTW.REDFT01, 1)
        else                        # n1, n2
            # iRFFT then iDST
            getproperty(state,qq)[2:end-1,:] = .5/(grid_y-2)*FFTW.r2r(
                real(
                    FFTW.irfft(
                        getproperty(statehat,qq), grid_x, 2)
                    )[2:end,:],
                FFTW.RODFT01, 1)
        end
    end
    return state
end


#=Drops last Coeff for DCT/DST
function dcsft(state::MJO_State, statehat::MJO_State_im,
    fields=fieldnames(MJO_State)[1:6], grid_x::Int=1440)
    for qq in fields
        if qq in fields[[1,3,5,6]]   # m1, m2, h1, h2 
            # DCT then RFFT
            getproperty(statehat,qq)[:,:] = FFTW.rfft(
                FFTW.r2r(getproperty(state,qq)[2:end-1,:], FFTW.REDFT10, 1)
                #=vcat(
                    FFTW.r2r(getproperty(state,qq)[2:end-1,:], FFTW.REDFT10, 1),
                    zeros(1,grid_x)
                    )=#, 2)
        else                        # n1, n2
            # DST then RFFT
            getproperty(statehat,qq)[:,:] = FFTW.rfft(
                vcat(
                    zeros(1,grid_x),
                    FFTW.r2r(getproperty(state,qq)[2:end-1,:], FFTW.RODFT10, 1)[1:end-1, :]
                    ), 2)
        end
    end
    return statehat
end

function idcsft(state::MJO_State, statehat::MJO_State_im,
    fields=fieldnames(MJO_State)[1:6], grid_x::Int=1440)
    for qq in fields
        if qq in fields[[1,3,5,6]]   # m1, m2, h1, h2 
            # iRFFT then iDCT
            getproperty(state,qq)[2:end-1,:] = .5/(grid_y-2)*FFTW.r2r(
                real(
                    FFTW.irfft(getproperty(statehat,qq), grid_x, 2)
                    ),#[1:end-1,:],
                FFTW.REDFT01, 1)
        else                        # n1, n2
            # iRFFT then iDST
            getproperty(state,qq)[2:end-1,:] = .5/(grid_y-2)*FFTW.r2r(
                vcat(
                    real(
                    FFTW.irfft(
                        getproperty(statehat,qq)
                        , grid_x, 2)
                    )[2:end,:],
                    zeros(1,grid_x)
                    ),
                FFTW.RODFT01, 1)
        end
    end
    return state
end=#
