
abstract type AbstractProblemParSetter end
struct DummyProblemParSetter <: AbstractProblemParSetter end # for testing error message

"""
    count_state(::AbstractProblemParSetter) 
    count_par(::AbstractProblemParSetter) 
    count_paropt(::AbstractProblemParSetter) 

Report the number of problem states, problem parameters and optimized parameters
respectively.    
"""
function count_state end, function count_par end, function count_paropt end

@deprecate count_states(pset::AbstractProblemParSetter) count_state(pset)

"""
    symbols_state(pset::AbstractProblemParSetter)
    symbols_par(pset::AbstractProblemParSetter)
    symbols_paropt(pset::AbstractProblemParSetter)

Report the names, i.e. symbols of problem states, problem parameters and 
optimized parameters respectively, i.e. the concatenation of components.
Similar to `ComponentArrays.label`, but inferred from Axis object.   
"""
function symbols_state end, function symbols_par end, function symbols_paropt end
# need to implement in concrete types

@deprecate statesyms(pset::AbstractProblemParSetter) symbols_state(pset)
@deprecate parsyms(pset::AbstractProblemParSetter) symbols_par(pset)
@deprecate paroptsyms(pset::AbstractProblemParSetter) symbols_paropt(pset)


"""
    axis_state(pset::AbstractProblemParSetter)
    axis_par(pset::AbstractProblemParSetter)
    axis_paropt(pset::AbstractProblemParSetter)

Report the Axis, i.e. nested component symbols of problem states, problem parameters and 
optimized parameters respectively.
Returns an AbstractAxis.
"""
function axis_state(pset::AbstractProblemParSetter) end,
function axis_par(pset::AbstractProblemParSetter) end,
function axis_paropt(pset::AbstractProblemParSetter) end
# need to implement in concrete types


"""
    get_paropt(pset::AbstractProblemParSetter, prob::ODEProblem; kwargs...)
    get_paropt(pset::AbstractProblemParSetter, u0, p)

    get_paropt_labeled(pset::AbstractProblemParSetter, prob::ODEProblem; kwargs...)
    get_paropt_labeled(pset::AbstractProblemParSetter, u0, p)

Extract the initial states and parameters corresponding to the positions
that are optimized.    
If both u0 and p are AbstractVectors, the result is a Vector, otherwise the result is a Tuple.

The labeled versions additionally call `label_paropt` (see [`label_state`](@ref)) 
on the return value.
"""
function get_paropt(pset::AbstractProblemParSetter, prob::ODEProblem; kwargs...)
    get_paropt(pset, prob.u0, prob.p; kwargs...)
end,
function get_paropt_labeled(pset::AbstractProblemParSetter, prob::ODEProblem; kwargs...)
    get_paropt_labeled(pset, prob.u0, prob.p; kwargs...)
end,
function get_paropt(pset::AbstractProblemParSetter, u0, p)
    getdata(get_paropt_labeled(pset, u0, p))

end
# need to implement in concrete types: get_paropt_labeled -> ComponentVector

"""
    label_state(pset::AbstractProblemParSetter, u::AbstractVector) 
    label_par(pset::AbstractProblemParSetter, par::AbstractVector) 
    label_paropt(pset::AbstractProblemParSetter, popt::AbstractVector) 

Produce a labeled version, i.e. a ComponentVector of initial states, parameters, or
optimized parameters respectively.
"""
function label_state(pset::AbstractProblemParSetter, u)
    attach_axis(u, axis_state(pset))
end,
function label_par(pset::AbstractProblemParSetter, p)
    attach_axis(p, axis_par(pset))
end,
function label_paropt(pset::AbstractProblemParSetter, popt)
    attach_axis(popt, axis_paropt(pset))
end

# TODO move to ComponentArrays.jl
# type piracy: https://github.com/jonniedie/ComponentArrays.jl/issues/141
#@inline CA.getdata(x::ComponentVector) = getfield(x, :data)
attach_axis(x::AbstractVector, ax::AbstractAxis) = ComponentArray(x, (ax,))
#attach_axis(x::ComponentVector, ax::AbstractAxis) = ComponentArray(getdata(x), (ax,))
attach_axis(x::ComponentVector, ax::AbstractAxis) =
    ComponentArray(getfield(x, :data), (ax,))
attach_x_axis(x::ComponentMatrix, ax::AbstractAxis) = ComponentArray(x, (ax, FlatAxis()))


# label_state(pset::AbstractProblemParSetter, u::SVector) = SLVector(label_state(pset, Tuple(u)))
# label_state(pset::AbstractProblemParSetter, u::NTuple) = NamedTuple{symbols_state(pset)}(u)

# label_par(pset::AbstractProblemParSetter, par::SVector) = SLVector(label_par(pset, Tuple(par)))
# label_par(pset::AbstractProblemParSetter, par::NTuple) = NamedTuple{symbols_par(pset)}(par)

# label_paropt(pset::AbstractProblemParSetter, popt::SVector) = SLVector(label_paropt(pset, Tuple(popt)))
# label_paropt(pset::AbstractProblemParSetter, popt::NTuple) = NamedTuple{symbols_paropt(pset)}(popt)

"""
    name_state(pset, u::AbstractVector) 
    name_par(pset, par::AbstractVector) 
    name_paropt(pset, popt::AbstractVector) 

    name_paropt(pset::AbstractProblemParSetter, prob::ODEProblem)    

Produce a `NamedVector` of given state, parameters, or optimized vars.
"""
function name_state(pset::AbstractProblemParSetter, state::AbstractVector)
    NamedArray(state, (collect(symbols_state(pset))::Vector{Symbol},))
end,
function name_par(pset::AbstractProblemParSetter, par::AbstractVector)
    NamedArray(par, (collect(symbols_par(pset))::Vector{Symbol},))
end,
function name_paropt(pset::AbstractProblemParSetter, paropt::AbstractVector)
    NamedArray(paropt, (collect(symbols_paropt(pset))::Vector{Symbol},))
end,
function name_paropt(pset::AbstractProblemParSetter, prob::ODEProblem; kwargs...)
    name_paropt(pset, get_paropt(pset, prob); kwargs...)
end

"""
    prob_new = update_statepar(ps::ProblemParSetter, popt, prob::ODEProblem) 
    u0new, pnew = update_statepar(ps::ProblemParSetter, popt, u0, p) 

Return an updated problem or updates states and parameters where
values corresponding to positions in `popt` hold the values of popt.
The type is changed to the promotion type of popt, to allow working with dual numbers.
"""
function update_statepar(pset::AbstractProblemParSetter, popt, prob::ODEProblem)
    u0, p = update_statepar(pset, popt, prob.u0, prob.p)
    remake(prob; u0, p)
end
