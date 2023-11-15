module MTKHelpers

using ModelingToolkit, OrdinaryDiffEq
using ModelingToolkit: AbstractODESystem
using SciMLBase: SciMLBase, AbstractODEProblem
using StaticArrays, LabelledArrays
using NamedArrays: NamedArrays
using ComponentArrays
using Distributions
using SymbolicUtils: SymbolicUtils
using InlineStrings
#using Infiltrator

export AbstractProblemParSetter,
    AbstractODEProblemParSetter,
    ODEProblemParSetter,
    count_state,
    count_par,
    count_paropt,
    axis_paropt,
    axis_par,
    axis_state,
    keys_paropt,
    keys_par,
    keys_state,
    symbols_paropt,
    symbols_state,
    symbols_par,
    update_statepar,
    get_paropt,
    get_paropt_labeled,
    label_paropt,
    label_par,
    label_state,
    name_paropt,
    name_par,
    name_state,
    get_u_map,
    get_p_map


export smoothstep

export series_sol!

export getlast

# extending 
import Base: merge, getindex
import ComponentArrays: getdata

export symbol_op, symbols_state, symbols_par, strip_namespace, embed_system, override_system
include("util.jl")

#export _get_index_axis, _set_index_axis!, attach_axis, _update_cv, _labels
include("util_componentarrays.jl")

include("abstractproblemparsetter.jl")
include("abstractodeproblemparsetter.jl")
include("odeproblemparsetter.jl")
include("smoothstep.jl")
include("solution.jl")

export AbstractProblemUpdater,
    AbstractProblemParGetter, ProblemUpdater, KeysProblemParGetter, NullProblemUpdater
export par_setter, par_setter
export get_ode_problemupdater
include("problem_updater.jl")

export get_system_symbol_dict, system_num_dict
export strip_deriv_num
include("util_nums.jl")

export fit_Dirichlet_std, fit_Dirichlet_mode, simplex_grid
include("prior_util.jl")

export series_sol!
include("makie_util.jl")

#export samplesystem_vec, indices_of_nums
include("example_systems.jl")


if !isdefined(Base, :get_extension)
    using Requires
end

@static if !isdefined(Base, :get_extension)
    function __init__()
        @require CairoMakie="13f3f980-e62b-5c42-98c6-ff1f3baf88f0" include("../ext/MTKHelpersMakieExt.jl")
    end
end

end # module
