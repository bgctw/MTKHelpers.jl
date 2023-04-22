module MTKHelpers

using ModelingToolkit, DifferentialEquations
using StaticArrays, LabelledArrays
using NamedArrays
using ComponentArrays
using Distributions
#using Infiltrator

export AbstractProblemParSetter, ProblemParSetter_sym, ProblemParSetter,
    count_state, count_par, count_paropt, 
    axis_paropt, axis_par, axis_state,
    symbols_paropt, symbols_state, symbols_par,
    update_statepar, get_paropt, get_paropt_labeled, label_paropt, label_par, label_state,
    name_paropt, name_par, name_state,
    get_u_map, get_p_map

export    
    _get_index_axis, _set_index_axis!, attach_axis, _update_cv, _labels

export smoothstep

export series_sol!

export getlast

# extending 
import Base: merge, getindex    
import ComponentArrays: getdata

if !isdefined(Base, :get_extension)
    using Requires
end

@static if !isdefined(Base, :get_extension)
    function __init__()
        @require CairoMakie = "13f3f980-e62b-5c42-98c6-ff1f3baf88f0" include("../ext/MTKHelpersMakieExt.jl")
    end
end

"""
    series_sol!(ax, sol::AbstractODESolution, vars; tspan=extrema(sol.t), labels=string.(vars), nt=120, kwargs...)

calls `CairoMakie.series` for a grid fo `n` points and interpolated values
from `sol`.
Currently works only with solutions created by a non-composite solver, e.g. `Tsit5`.
"""
function series_sol! end
export series_sol!
  
export symbol, symbols_state, symbols_par, strip_namespace, embed_system, override_system
include("util.jl")

include("util_componentarrays.jl")
include("abstractproblemparsetter.jl")
include("problemparsetter_sym.jl")
include("problemparsetter.jl")
include("smoothstep.jl")
include("solution.jl")

export AbstractProblemUpdater, AbstractProblemParGetter, ProblemUpdater, KeysProblemParGetter, NullProblemUpdater
include("problem_updater.jl")

export get_system_symbol_dict, system_num_dict
include("util_nums.jl")

export fit_Dirichlet_std, fit_Dirichlet_mode, simplex_grid
include("prior_util.jl")

end
