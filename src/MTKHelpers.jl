module MTKHelpers

using ModelingToolkit, DifferentialEquations
using StaticArrays, LabelledArrays
using Requires: @require 
using NamedArrays
#using Infiltrator

export symbol, symbols_state, symbols_par, strip_namespace, embed_system, cm2inch

export ProblemParSetter, count_state, count_par, count_paropt, 
    symbols_paropt, symbols_state, symbols_par,
    update_statepar, get_paropt, get_paropt_labeled, label_paropt, label_par, label_state,
    name_paropt, name_par, name_state

export smoothstep

export series_sol!

export getlast

# extending 
import Base: merge    

function __init__()
    @require CairoMakie = "13f3f980-e62b-5c42-98c6-ff1f3baf88f0" include("requires_cairommakie.jl")
  end
  
  
include("util.jl")
include("problemparsetter.jl")
include("smoothstep.jl")
include("solution.jl")



end
