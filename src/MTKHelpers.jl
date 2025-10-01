module MTKHelpers

using ModelingToolkit, OrdinaryDiffEq
using ModelingToolkit: t_nounits as t, D_nounits as D
using SciMLBase: SciMLBase, AbstractODEProblem
using StaticArrays, LabelledArrays
using NamedArrays: NamedArrays
using ComponentArrays
using Distributions
using SymbolicUtils: SymbolicUtils
using InlineStrings
using Chain
using LoggingExtras: LoggingExtras
using Logging: Logging
using SciMLStructures
using SymbolicIndexingInterface: setp, SymbolicIndexingInterface as SII
#using Infiltrator

export AbstractProblemParSetter,
    AbstractODEProblemParSetter,
    ODEProblemParSetter,
    ODEProblemParSetterConcrete,
    count_state,
    count_par,
    count_paropt,
    axis_paropt,
    axis_paropt_scalar,
    axis_paropt_flat1,
    axis_par,
    axis_state,
    axis_state_scalar,
    keys_paropt,
    keys_par,
    keys_state,
    symbols_paropt,
    symbols_state,
    symbols_par,
    #update_statepar,
    get_paropt,
    get_paropt_labeled,
    label_paropt,
    label_paropt_flat1,
    get_par,
    get_par_labeled,
    get_state,
    get_state_labeled,
    label_par,
    label_state,
    name_paropt,
    name_par,
    name_state,
    #get_u_map,
    #get_p_map,
    get_concrete

# extending 
import Base: merge, getindex
import ComponentArrays: getdata

export symbol_op, embed_system, override_system
include("util.jl")

#export _get_index_axis, _set_index_axis!, attach_axis, _update_cv, _labels
export flatten1
export vcat_statesfirst
#export map_keys
include("util_componentarrays.jl")

include("abstractproblemparsetter.jl")

export get_system
include("abstractodeproblemparsetter.jl")

export NullODEProblemParSetter, NullODEProblemParSetterConcrete
include("nullodeproblemparsetter.jl")

include("odeproblemparsetterconcrete.jl")
include("odeproblemparsetter.jl")


export smoothstep
include("smoothstep.jl")

export getlast
include("solution.jl")

export AbstractProblemUpdater, ProblemUpdater, ProblemUpdaterConcrete, NullProblemUpdater,
    AbstractProblemParGetter, KeysProblemParGetter
export par_setter, par_setter
export get_ode_problemupdater
include("problemupdater.jl")

export get_system_symbol_dict, system_num_dict
export strip_deriv_num
export base_num
#export get_base_num_dict
include("util_nums.jl")

export fit_Dirichlet_std, fit_Dirichlet_mode, simplex_grid
include("prior_util.jl")

export grid_exp, Dz_exp, Iz_exp, Dz_lin, Iz_lin
export get_1d_state_pos, get_1d_grid, get_discrete_space
include("util_pde.jl")

export series_sol!
include("makie_util.jl")

#export samplesystem_vec, indices_of_nums
include("example_systems.jl")

export read_csv_cv, write_csv_cv
include("cvwriter.jl")

if !isdefined(Base, :get_extension)
    using Requires
end

@static if !isdefined(Base, :get_extension)
    function __init__()
        @require CairoMakie="13f3f980-e62b-5c42-98c6-ff1f3baf88f0" include("../ext/MTKHelpersMakieExt.jl")
        @require MethodOfLines="94925ecb-adb7-4558-8ed8-f975c56a0bf4" include("../ext/MTKHelpersMethodOfLinesExt.jl")
        @require CSV = "336ed68f-0bac-5ca0-87d4-7b16caf5d00b"
         include("../ext/MTKHelpersCSVExt.jl")
    end
end

end # module
