function smoother(q::Array{Float64,2}, stencil::Array{Int,2}; 
    periodic::Bool=true)
    m, n = size(stencil);
    mm = floor(Int,m/2)
    nn = floor(Int,n/2)
    grid_y, grid_x = size(q)
    smoothed = Array{Float64,2}(undef, grid_y, grid_x)
    if periodic==false
        for i = 1 : grid_y
            for j = 1 : grid_x
                top_s, bottom_s, left_s, right_s = 1,m,1,n;
                top_q, bottom_q, left_q, right_q = i-mm,i+mm,j-nn,j+nn;
                if i<=mm
                    top_s=mm-i+2;           top_q=1;
                elseif i>grid_y-mm
                    bottom_s=mm+grid_y-i+1; bottom_q=grid_y;
                end
                if j<=nn
                    left_s=nn-j+2;          left_q=1;
                elseif j>grid_x-nn
                    right_s=nn+grid_x-j+1;  right_q=grid_x;
                end
                smoothed[i,j] = 1/sum(stencil[top_s:bottom_s, left_s:right_s])*(
                    sum(q[top_q:bottom_q, left_q:right_q]
                        .*stencil[top_s:bottom_s, left_s:right_s])
                    )
            end
        end
        return smoothed
    else
        q = hcat(q[:,end-nn+1:end], q, q[:, 1:nn])
        for i = 1 : grid_y
            for j = nn+1 : grid_x-nn
                top_s, bottom_s, left_s, right_s = 1,m,1,n;
                top_q, bottom_q, left_q, right_q = i-mm,i+mm,j-nn,j+nn;
                if i<=mm
                    top_s=mm-i+2;           top_q=1
                elseif i>grid_y-mm
                    bottom_s=mm+grid_y-i+1; bottom_q=grid_y
                end
                smoothed[i,j] = 1/sum(stencil[top_s:bottom_s, left_s:right_s])*(
                    sum(q[top_q:bottom_q, left_q:right_q]
                        .*stencil[top_s:bottom_s, left_s:right_s])
                    )
            end
        end
        return smoothed
    end
end