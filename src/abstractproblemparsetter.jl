"""
During an optimization, one does not want to recreate the problem
from a symbolic system, but only update the problem.
This can be difficult, because the parameters to update may be spread
across initial state and other parameters. 
Further, the order of parameters after simplifying a system is not fixed.

A AbstractProblemUpdater helps with 
- [`remake`](@ref): translate the set of parameters -> an updated problem 
- [`get_paropt`](@ref): problem -> extract/approximate subset of parameters to optimize 

The structure of optimized parameter Vector is described by an Axis object 
of ComponentArrays.jl. And several functions are defined to work with it.
Specifically, the ComponentVector it employs a classification ([`classes_paropt`](@ref)),
,e.g. :`state` and `:par` for ODEProblems,
below which, ComponentVectors of actual parameters are listed (`keys_paropt`](@ref)).
Further functions extract information about the ComponentVector:
[`axis_paropt`](@ref), , [`count_paropt`](@ref), 
[`symbols_paropt`](@ref)
or attach information to a plain vector for convenient access or display:
[`label_paropt`](@ref), [`name_paropt`](@ref). 
"""
abstract type AbstractProblemParSetter end

"""
    remake(prob::AbstractSciMLProblem, popt, ps::AbstractProblemParSetter) 

Return an updated problem given the parameters.
Subtypes need to implement method `remake_pset(prob, popt, pset)`
"""
function SciMLBase.remake(prob::SciMLBase.AbstractSciMLProblem, popt, pset::AbstractProblemParSetter) 
    remake_pset(prob, popt, pset)
end

"""
    get_paropt(pset::AbstractProblemParSetter, prob::SciMLBase.AbstractSciMLProblem; kwargs...)
    get_paropt_labeled(pset::AbstractProblemParSetter, prob::SciMLBase.AbstractSciMLProblem; kwargs...)

Extract optimized parameters from the Problem.
The labeled versions additionally calls [`label_paropt`](@ref)
on the return value.
"""
function get_paropt(pset::AbstractProblemParSetter, prob::SciMLBase.AbstractSciMLProblem; kwargs...) end,
function get_paropt_labeled(pset::AbstractProblemParSetter, prob::SciMLBase.AbstractSciMLProblem; kwargs...) end


"""
    axis_paropt(pset::AbstractProblemParSetter)

Report the Axis of a CompoenentVector of parameters.
"""
function axis_paropt(::AbstractProblemParSetter) end

"""
    function classes_paropt(pset::AbstractProblemParSetter) 

Get the classes (as NTuple{Symbol}) which the AbstractProblemParSetter 
supports and requires in paropt.        
"""
function classes_paropt(::AbstractProblemParSetter) end
# need to implement in concrete types

"""
    count_paropt(::AbstractProblemParSetter) 

Report length of the optimized parameters vector.    
This generally is different from the length of keys, because each key
can describe a array.
"""
function count_paropt(pset::AbstractProblemParSetter) 
    axis_length(axis_paropt(pset))::Int
end
# TODO change to length(ax) when this becomes available in ComponentArrays

"""
    keys_paropt(::AbstractProblemParSetter) 

Report the keys of paropt below the classification level.
"""
function keys_paropt(ps::AbstractProblemParSetter) 
    ax = axis_paropt(ps)
    # for each top-key access the subaxis and apply keys
    gen = (getproperty(CA.indexmap(ax), k) |> x -> keys(x) for k in classes_paropt(ps))
    tuplejoin(gen...)
end

"""
    symbols_paropt(pset::AbstractProblemParSetter)

Report the names, i.e. symbols of optimiz parameters respectively, 
i.e. the concatenation of components.
Similar to `ComponentArrays.label`, but inferred from Axis object.  
Returns a Vector of length [`count_paropt`](@ref)
"""
function symbols_paropt(pset::AbstractProblemParSetter)
    ax = axis_paropt(pset)
    # for each key access the subaxis and apply _ax_symbols_tuple
    gen = (getproperty(CA.indexmap(ax), k) |> x -> _ax_symbols_tuple(x)
           for k in classes_paropt(pset))
    # concatenate the generator of tuples
    tuplejoin(gen...)
end
#symbols_paropt(pset::AbstractProblemParSetter) = _ax_symbols_tuple(axis_paropt(pset))
# concatenate the symbols of subaxes

@deprecate paroptsyms(pset::AbstractProblemParSetter) symbols_paropt(pset)

"""
    label_paropt(pset::AbstractProblemParSetter, popt::AbstractVector) 

Produce a labeled version, i.e. a ComponentVector of optimized parameters.
"""
function label_paropt(pset::AbstractProblemParSetter, popt)
    attach_axis(popt, axis_paropt(pset))
end


"""
    name_paropt(pset, popt::AbstractVector) 
    name_paropt(pset, prob::AbstractSciMLProblem)    

Produce a `NamedVector` of given state optimized parameters vector.
Similar to [`label_paropt`](@ref), but may print nicer.
The second form calls [`get_paropt`](@ref) on the Problem.
"""
function name_paropt(pset::AbstractProblemParSetter, paropt::AbstractVector)
    NamedArrays.NamedArray(paropt, (collect(symbols_paropt(pset))::Vector{Symbol},))
end,
function name_paropt(pset::AbstractProblemParSetter, prob::SciMLBase.AbstractSciMLProblem; kwargs...)
    name_paropt(pset, get_paropt(pset, prob); kwargs...)
end

"""
    get_concrete(pset::AbstractProblemParSetter)
    get_concrete(pu::AbstractParameterUpdater)

Return a concrete-type-version of an ProblemParSetter or ProblemUpdater.
"""
function get_concrete(pset::AbstractProblemParSetter)
    !isconcrete(pset) && @warn "no concrete type implemented for ProblemSetter of type $(typeof(pset))."
    pset
end













