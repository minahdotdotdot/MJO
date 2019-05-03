include("mjo_a.jl");
include("mjo_imex2.jl");
include("smooth_data.jl")
ENV["MPLBACKEND"]="Agg"
using PyPlot, Printf

function savecontourmaps(evol::Array{MJO_State,1}, str::String; 
    draw::Symbol=:contourf,
    aspect::Bool=true,
    loc::String="../movies/")
    for f in fieldnames(MJO_State)
        evolfield = 1
        cm = "PuOr"
        if f == :m1 || f ==:n1
            evolfield = elemdiv(getproperty(evol, f), getproperty(evol, :h1))
        elseif f == :m2 || f ==:n2
            evolfield = elemdiv(getproperty(evol, f), getproperty(evol, :h2))
        else
            evolfield = getproperty(evol, f)
            cm = "BuGn"
        end
        minval = minimum(evolfield);
        maxval = maximum(evolfield);
        for j = 1 : length(evolfield)
            fig, ax = subplots(); 
            fig[:set_size_inches](13,2); 
            if aspect==true
                ax[:set_aspect]("equal");
            end
            fig[:colorbar](
                ax[draw](
                    params.lon, 
                    params.lat[2:end-1],  
                    evolfield[j],
                    vmin=minval,
                    vmax=maxval,
                    cmap=cm
                )
            );
            savefig(
                loc*string(f)*"/"*str*string(j),
                pad_inches=.10, 
                bbox_inches="tight"
            )
            close(fig)
        end
    end
end

@inline function savecontourf(state::MJO_State, ii::String; loc::String="../movies/")
    for f in fieldnames(MJO_State)
        evolfield = 1
        cm = "PuOr"
        if f == :m1 || f ==:n1
            evolfield = getproperty(state,f)[2:end-1, :]./(1 .+getproperty(state,:h1)[2:end-1, :]);
        elseif f ==:m2 || f ==:n2
            evolfield = getproperty(state,f)[2:end-1, :]./(1 .+getproperty(state,:h2)[2:end-1, :]);
        else
            evolfield = getproperty(state,f)[2:end-1, :];
            if f ==:q
                cm = "BuGn"
            end
        end
        fig, ax = subplots(); 
        getproperty(fig, :set_size_inches)((13,2))
        getproperty(ax, :set_aspect)("equal")
        #fig[:set_size_inches](13,2); 
        #ax[:set_aspect]("equal");
        fig.colorbar(
                ax.contourf(
        #fig[:colorbar](
        #    ax[:contourf](
                params.lon,
                params.lat[2:end-1,:],
                evolfield,
                cmap=cm
            )
        );
        savefig(
            loc*string(f)*"/"*ii,
            pad_inches=.10, 
            bbox_inches="tight"
        )
        close(fig)
    end
end

@inline function saveimshow(state::MJO_State, ii::String; loc::String="../movies/", params::MJO_params, H1::Float64=1.0)
    H2 = 2.0 - H1
    for f in fieldnames(MJO_State)
        evolfield = 1
        cm = "PuOr"
        if f == :m1 || f ==:n1
            evolfield = getproperty(state,f)[2:end-1, :]./(H1 .+getproperty(state,:h1)[2:end-1, :]);
        elseif f ==:m2 || f ==:n2
            evolfield = getproperty(state,f)[2:end-1, :]./(H2 .+getproperty(state,:h2)[2:end-1, :]);
        elseif f==:q
            cm = "gist_ncar"#"BuGn"
            evolfield = P(params.LL, params.UU, params.QQ, params.B, params.Qs, getproperty(state,f)[2:end-1, :],
                params.T_Q, params.PP);
        else
            evolfield = getproperty(state,f)[2:end-1, :];
        end
        fig, ax = subplots(); 
        getproperty(fig, :set_size_inches)((10,6)) #(13,2)
        getproperty(ax, :set_aspect)("equal")
        #fig[:set_size_inches](13,2); 
        #ax[:set_aspect]("equal");
        fig.colorbar(
        	ax.imshow(
        #fig[:colorbar](
        #    ax[:imshow](
                evolfield,
                cmap=cm,
                extent=(params.lon_range[1], params.lon_range[2], params.lat_range[1], params.lat_range[2])#(0, 360, -20, 20)
            )
        );
        savefig(
            loc*string(f)*"/"*ii,
            pad_inches=.10, 
            bbox_inches="tight"
        )
        close(fig)
    end
end





# Should have net north momentum = 0 (integral/sum) 
# & zero at boundaries
# check periodic boundary and make sure initial conditions are periodic in x
# zero momenta and constant heightfields
# random moisture field (random, or random then smooth it. )
# show velocity instead of momenta DONE
# see how long it takes for interesting moisture dynamics
    # interesting == breaking code
    # or dynamics 
    # 20 days ish
# github
# future: RK4 (implicit maybe)

# thoughts: make P be a function of concentration (instead of total water)
# i.e. Q <- Q/h1, Qs a different appropriate parameter 
# try getting rid of div terms system of 3 equations




