module MTKHelpersMakieExt

function __init__()
    @info "MTKHelpers: loading MTKHelpersMakieExt"
end

isdefined(Base, :get_extension) ? (using CairoMakie) : (using ..CairoMakie)
using DiffEqBase: AbstractODESolution
using RecursiveArrayTools: VectorOfArray
import MTKHelpers # provide method for series_sol!

# pdf_figure moved to TwPrototypes

function MTKHelpers.series_sol!(ax,
        sol::AbstractODESolution,
        vars;
        tspan = extrema(sol.t),
        labels = string.(vars),
        nt = 120,
        kwargs...)
    #ts = first(tspan) .<= sol.t .<= last(tspan)
    #series!(ax, sol.t[ts], transpose(VectorOfArray(sol[vars])[ts,:]); labels, kwargs...)
    tsol = sol.t[first(tspan) .<= sol.t .<= last(tspan)]
    ts = range(first(tspan), last(tspan), length = nt)
    ts2 = sort!(vcat(tsol, ts))
    a = VectorOfArray(sol(ts2, idxs = vars).u)
    @show a
    series!(ax, ts2, a; labels, kwargs...)
end

end
