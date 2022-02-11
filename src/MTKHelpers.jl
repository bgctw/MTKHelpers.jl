module MTKHelpers

using ModelingToolkit, DifferentialEquations
using StaticArrays, LabelledArrays
using Requires: @require 


export symbol, statesyms, parsyms, strip_namespace, embed_system

export ProblemParSetter, count_states, count_par, count_paropt, 
    paroptsyms, statesyms, parsyms,
    update_statepar, get_paropt, get_paropt_labeled, label_paropt, label_par, label_state

export smoothstep

export pdf_figure, series_sol!

# extending 
import Base: merge    

function __init__()
    @require CairoMakie = "13f3f980-e62b-5c42-98c6-ff1f3baf88f0" include("requires_cairommakie.jl")
  end
  
  
include("util.jl")
include("problemparsetter.jl")
include("smoothstep.jl")



end
