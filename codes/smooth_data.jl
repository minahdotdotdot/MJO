function genStencil(m::Int, n::Int)
    # generates a gaussian kernel 
    if isodd.([m,n])==[true, true]
        x , y = 0, 0;
        stencil = Array{Float64, 2}(undef, m, n);
        if m > n
            x = range(-1, stop=1, length=n)
            y = vcat(
                range(0, step=-step(x), length=ceil(Int,m/2))[end:-1:2], 
                range(0, step=step(x), length=ceil(Int,m/2))
                )
        elseif m < n
            y = range(-1, stop=1, length=m)
            x = vcat(
                range(0, step=-step(y), length=ceil(Int,n/2))[end:-1:2], 
                range(0, step=step(y), length=ceil(Int,n/2))
                )
        else
            x = range(-1, stop=1, length=n)
            y = copy(x)
        end
        for j = 1 : m
            for i = 1 : n
                stencil[j,i] = exp(-x[i]^2-y[j]^2)
            end
        end
        return stencil/norm(stencil,2)
    end
end

function quadratic(a, b, c)
    discr = b^2 - 4*a*c
    if discr >= 0
        return (-b + sqrt(discr))/(2a), (-b - sqrt(discr))/(2a)
    else error("Only complex roots")
    end
end

function cov(;a::Float64=0.0, b::Float64=0.0)
    stencil = ones(3,3)
    if a ==0 && b==0
        error("Need to provide at least one var (a or b, 1>a>b)")
    elseif b ==0
        b1, b2 = quadratic(3, 2+8*a,2(a+a^2)-.25)
        if b1>0 || b2>0
            b = max(b1, b2)
        else
            error("no positive b values.")
        end
    elseif a ==0
        a1, a2 = quadratic(2, 2+8*b, b*(2+3*b)-.25)
        if a1>0 || a2>0
            a = max(a1, a2)
        else
            error("no positive a values.")
        end
    end
    for i = 1 : 3
        for j = 1 : 3
            if abs(i-j)==1
                stencil[i,j] = a
            else
                stencil[i,j] = b
            end
        end
    end
    stencil[2,2]=1.0
    C = [1+4*(a^2+b^2), 2*a*(1+2*b), 2*(b+a^2),a^2+2*b^2, 2*a*b, 2*b^2]/(1+4*(a^2+b^2))
    return stencil/sqrt(1+4*a^2+4*b^2)#, C
end 

function smoother(q::Array{T,2}; 
    stencil::Array{T,2}=cov(a=.11),
    periodic::Bool=true) where T<:Real
    m, n = size(stencil);
    mm = floor(Int,m/2)
    nn = floor(Int,n/2)
    grid_y, grid_x = size(q)
    smoothed = Array{Float64,2}(undef, grid_y, grid_x)
    if periodic==false
        for j = 1 : grid_y
            for i = 1 : grid_x
                top_s, bottom_s, left_s, right_s = 1,m,1,n;
                top_q, bottom_q, left_q, right_q = j-mm,j+mm,i-nn,i+nn;
                if j<=mm
                    top_s=mm-j+2;           top_q=1;
                elseif j>grid_y-mm
                    bottom_s=mm+grid_y-j+1; bottom_q=grid_y;
                end
                if i<=nn
                    left_s=nn-i+2;          left_q=1;
                elseif i>grid_x-nn
                    right_s=nn+grid_x-i+1;  right_q=grid_x;
                end
                smoothed[j,i] = 1/sum(stencil[top_s:bottom_s, left_s:right_s])*(
                    sum(q[top_q:bottom_q, left_q:right_q]
                        .*stencil[top_s:bottom_s, left_s:right_s])
                    )
            end
        end
        return smoothed
    else
        q = hcat(q[:,end-nn+1:end], q, q[:, 1:nn])
        for j = 1 : grid_y
            for i = nn+1 : grid_x+nn
                top_s, bottom_s, left_s, right_s = 1,m,1,n;
                top_q, bottom_q, left_q, right_q = j-mm,j+mm,i-nn,i+nn;
                if j<=mm
                    top_s=mm-j+2;           top_q=1
                elseif j>grid_y-mm
                    bottom_s=mm+grid_y-j+1; bottom_q=grid_y
                end
                smoothed[j,i-nn] = 1/sum(stencil[top_s:bottom_s, left_s:right_s])*(
                    sum(q[top_q:bottom_q, left_q:right_q]
                        .*stencil[top_s:bottom_s, left_s:right_s])
                    )
            end
        end
        return smoothed
    end
end

using GaussianRandomFields, PyPlot
function genRandField()
    #mat = Matern(0.5,2.0)
    #covar = CovarianceFunction(2,mat)
    grf = GaussianRandomField(
        CovarianceFunction(2,Matern(0.5,2.0)),
        CirculantEmbedding(),
        1:1:162, 
        1:1:1440,
        minpadding=250
        )
    return sample(grf)
end
