@info "MTKHelpers: loading CairoMakie utils"
using .CairoMakie # syntax by Requires.jl otherwise warning
using DiffEqBase: AbstractODESolution

"""
    pdf_figure = (size_inches = cm2inch.((8.3,8.3/1.618)); fontsize=9, kwargs...)

Creates a figure and Axis with specified resolution and fontsize 
for given figure size. See also [`cm2inch`](@ref) and `save`.
"""    
pdf_figure = (size_inches = cm2inch.((8.3,8.3/1.618)); fontsize=9, kwargs...) -> begin
    resolution = 72 .* size_inches # size_pt
    fig = Figure(;resolution, fontsize)
    fig, Axis(fig[1,1]; kwargs...)
end

"""
    series_sol!(ax, sol::AbstractODESolution, vars; tspan=extrema(sol.t), labels=string.(vars), nt=120, kwargs...)

calls `CairoMakie.series` for a grid fo `n` points and interpolated values
from `sol`.
Currently works only with solutions created by a non-composite solver, e.g. `Tsit5`.
"""
function series_sol!(ax, sol::AbstractODESolution, vars; tspan=extrema(sol.t), labels=string.(vars), nt=120, kwargs...)
    #ts = first(tspan) .<= sol.t .<= last(tspan)
    #series!(ax, sol.t[ts], transpose(VectorOfArray(sol[vars])[ts,:]); labels, kwargs...)
    tsol = sol.t[first(tspan) .<= sol.t .<= last(tspan)]
    ts = range(first(tspan), last(tspan), length=nt)
    ts2 = sort!(vcat(tsol,ts))
    series!(ax, ts2, VectorOfArray(sol(ts2, idxs=vars).u); labels, kwargs...)
end
