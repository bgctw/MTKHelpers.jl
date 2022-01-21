module MTKHelpers

using ModelingToolkit, DifferentialEquations
using StaticArrays, LabelledArrays

export symbol, statesyms, parsyms, strip_namespace, embed_system

export ProblemParSetter, count_states, count_par, count_paropt, 
    paroptsyms, statesyms, parsyms,
    update_statepar, get_paropt, label_paropt, label_par, label_state

# extending 
import Base: merge    


include("util.jl")
include("problemparsetter.jl")

end
