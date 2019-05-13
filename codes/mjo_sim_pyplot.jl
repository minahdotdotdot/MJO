include("mjo_a.jl");
include("mjo_imex2.jl");
include("smooth_data.jl")
ENV["MPLBACKEND"]="Agg"
using PyPlot, Printf

@inline function savecontourf(state::MJO_State, ii::String; loc::String="../movies/", params::MJO_params, H1::Float64=1.0)
    H2 = 2.0 - H1
    for f in fieldnames(MJO_State)
        evolfield = 1
        cm = "bwr"; #"PuOr" 
        vmin = nothing; vmax = nothing; #fix colorbar for q/P field.
        if f == :m1 || f ==:n1
            evolfield = getproperty(state,f)[2:end-1, :]./(H1 .+getproperty(state,:h1)[2:end-1, :]);
        elseif f ==:m2 || f ==:n2
            evolfield = getproperty(state,f)[2:end-1, :]./(H2 .+getproperty(state,:h2)[2:end-1, :]);
        elseif f==:q
            cm = "gist_ncar"
            evolfield = P(params.LL, params.UU, params.QQ, params.B, params.Qs, getproperty(state,f)[2:end-1, :],
                params.T_Q, params.PP);
            vmin=0.0;
            vmax=20.0;
        else
            evolfield = getproperty(state,f)[2:end-1, :];
        end
        fig, ax = subplots(); 
        getproperty(fig, :set_size_inches)((10,6)) #(13,2)
        getproperty(ax, :set_aspect)("equal")
        fig.colorbar(
            ax.contourf(
                evolfield,
                cmap=cm,
                extent=(params.lon_range[1], params.lon_range[2], params.lat_range[1], params.lat_range[2]),
                vmin=vmin, vmax=vmax
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
        cm = "bwr"; #"PuOr" 
        vmin = nothing; vmax = nothing; #fix colorbar for q/P field.
        if f == :m1 || f ==:n1
            evolfield = getproperty(state,f)[2:end-1, :]./(H1 .+getproperty(state,:h1)[2:end-1, :]);
        elseif f ==:m2 || f ==:n2
            evolfield = getproperty(state,f)[2:end-1, :]./(H2 .+getproperty(state,:h2)[2:end-1, :]);
        elseif f==:q
            cm = "gist_ncar"
            evolfield = P(params.LL, params.UU, params.QQ, params.B, params.Qs, getproperty(state,f)[2:end-1, :],
                params.T_Q, params.PP);
            vmin=0.0;
            vmax=20.0;
        else
            evolfield = getproperty(state,f)[2:end-1, :];
        end
        fig, ax = subplots(); 
        getproperty(fig, :set_size_inches)((10,6)) #(13,2)
        getproperty(ax, :set_aspect)("equal")
        fig.colorbar(
        	ax.imshow(
                evolfield,
                cmap=cm,
                extent=(params.lon_range[1], params.lon_range[2], params.lat_range[1], params.lat_range[2]),
                vmin=vmin, vmax=vmax
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