"""
    series_sol!(ax, sol::AbstractODESolution, vars; tspan=extrema(sol.t), labels=string.(vars), nt=120, kwargs...)

calls `CairoMakie.series` for a grid fo `n` points and interpolated values
from `sol`.
`vars` is passed to the solution object and can contain observables.
Currently works only with solutions created by a non-composite solver, e.g. `Tsit5`.
"""
function series_sol! end
