include("mjo_a.jl");
include("smooth_data.jl")
using PyPlot, Printf

function savecontourmaps(evol::Array{MJO_State,1}, str::String; 
    draw::Symbol=:contourf,
    aspect::Bool=true)
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
                "../movies/"*string(f)*"/"*str*string(j),
                pad_inches=.10, 
                bbox_inches="tight"
            )
            close(fig)
        end
    end
end

@inline function savecontourf(state::MJO_State, ii::String)
    for f in fieldnames(MJO_State)
        evolfield = 1
        cm = "PuOr"
        if f == :m1 || f ==:n1
            evolfield = getproperty(state,f)[2:end-1, :]./getproperty(state,:h1)[2:end-1, :];
        elseif f ==:m2 || f ==:n2
            evolfield = getproperty(state,f)[2:end-1, :]./getproperty(state,:h2)[2:end-1, :];
        else
            evolfield = getproperty(state,f)[2:end-1, :];
            if f ==:q
                cm = "BuGn"
            end
        end
        fig, ax = subplots(); 
        fig[:set_size_inches](13,2); 
        ax[:set_aspect]("equal");
        fig[:colorbar](
            ax[:contourf](
                params.lon,
                params.lat[2:end-1,:],
                evolfield,
                cmap=cm
            )
        );
        savefig(
            "../movies/"*string(f)*"/"*ii,
            pad_inches=.10, 
            bbox_inches="tight"
        )
        close(fig)
    end
end

@inline function saveimshow(state::MJO_State, ii::String)
    for f in fieldnames(MJO_State)
        evolfield = 1
        cm = "PuOr"
        if f == :m1 || f ==:n1
            evolfield = getproperty(state,f)[2:end-1, :]./getproperty(state,:h1)[2:end-1, :];
        elseif f ==:m2 || f ==:n2
            evolfield = getproperty(state,f)[2:end-1, :]./getproperty(state,:h2)[2:end-1, :];
        else
            evolfield = getproperty(state,f)[2:end-1, :];
            if f ==:q
                cm = "BuGn"
            end
        end
        fig, ax = subplots(); 
        fig[:set_size_inches](13,2); 
        ax[:set_aspect]("equal");
        fig[:colorbar](
            ax[:imshow](
                evolfield,
                cmap=cm,
                extent=(0, 360, -20, 20)
            )
        );
        savefig(
            "../movies/"*string(f)*"/"*ii,
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




