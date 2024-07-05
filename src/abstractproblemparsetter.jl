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
function SciMLBase.remake(prob::SciMLBase.AbstractSciMLProblem,
        popt, pset::AbstractProblemParSetter)
    remake_pset(prob, popt, pset)
end

"""
    get_paropt(pset::AbstractProblemParSetter, prob::SciMLBase.AbstractSciMLProblem; kwargs...)
    get_paropt_labeled(pset::AbstractProblemParSetter, prob::SciMLBase.AbstractSciMLProblem; kwargs...)

Extract optimized parameters from the Problem.
The labeled versions additionally calls [`label_paropt`](@ref)
on the return value.
"""
function get_paropt(pset::AbstractProblemParSetter, prob::SciMLBase.AbstractSciMLProblem;
        kwargs...) end,
function get_paropt_labeled(pset::AbstractProblemParSetter,
        prob::SciMLBase.AbstractSciMLProblem;
        kwargs...) end


"""
    get_par(pset::AbstractProblemParSetter, prob::SciMLBase.AbstractSciMLProblem; kwargs...)
    get_par_labeled(pset::AbstractProblemParSetter, prob::SciMLBase.AbstractSciMLProblem; kwargs...)

Extract parameters in the order of parameters(sys) from the Problem.
The labeled versions additionally calls [`label_par`](@ref)
on the return value.
"""
function get_par(pset::AbstractProblemParSetter, prob::SciMLBase.AbstractSciMLProblem;
        kwargs...) end,
function get_par_labeled(pset::AbstractProblemParSetter,
        prob::SciMLBase.AbstractSciMLProblem;
        kwargs...) 
    p = get_par(pset, prob; kwargs...)
    label_par(pset, p)
end

"""
    axis_paropt(pset::AbstractProblemParSetter)
    axis_paropt_scalar(pset::AbstractProblemParSetter)
    axis_paropt_flat1(pset::AbstractProblemParSetter)

Report the Axis of a CompoenentVector of parameters.
The second version has a scalarized entry for state for each subvector of state. The third version provides an axis corresponding to 
`flatten1(paropt)`.
"""
function axis_paropt(::AbstractProblemParSetter) end
function axis_paropt_scalar(::AbstractProblemParSetter) end
function axis_paropt_flat1(::AbstractProblemParSetter) end

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
function keys_paropt(pset::AbstractProblemParSetter)
    ax = axis_paropt(pset)
    # for each top-key access the subaxis and apply keys
    gen = (getproperty(CA.indexmap(ax), k) |> x -> keys(x) for k in classes_paropt(pset))
    tuplejoin(gen...)
end

"""
    symbols_paropt(pset::AbstractProblemParSetter)

Report the names, i.e. symbols of optimized parameters respectively, 
i.e. the concatenation of components.
Similar to `ComponentArrays.label`, but inferred from Axis object.  
Returns a Vector of length [`count_paropt`](@ref)
"""
function symbols_paropt(pset::AbstractProblemParSetter)
    ax = axis_paropt_scalar(pset)
    # for each key access the subaxis and apply _ax_symbols_tuple
    #k = first(classes_paropt(pset))
    #x = getproperty(CA.indexmap(ax), k)
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
    label_paropt_flat1(pset::AbstractProblemParSetter, popt::AbstractVector) 

Produce a labeled version, i.e. a ComponentVector of optimized parameters.
The second version omits the highest level of labels, e.g. state and par in 
`ODEProblemParSetter`.
"""
function label_paropt(pset::AbstractProblemParSetter, popt)
    attach_axis(popt, axis_paropt(pset))
end,
function label_paropt_flat1(pset::AbstractProblemParSetter, popt)
    # use stored flat axis to avoid repeated reduce(vcat(ComponentVector 
    # with flatten1
    ax = axis_paropt_flat1(pset)
    ax isa FlatAxis && 
    @warn("Called label_paropt_flat1 with a ProblemParSetter that has no flat version.")
    attach_axis(popt, ax)
    #flatten1(label_paropt(pset, popt))
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
function name_paropt(pset::AbstractProblemParSetter, prob::SciMLBase.AbstractSciMLProblem;
        kwargs...)
    name_paropt(pset, get_paropt(pset, prob); kwargs...)
end

"""
    get_concrete(pset::AbstractProblemParSetter)
    get_concrete(pu::AbstractParameterUpdater)

Return a concrete-type-version of an ProblemParSetter or ProblemUpdater.
"""
function get_concrete(pset::AbstractProblemParSetter)
    !isconcrete(pset) &&
        @warn "no concrete type implemented for ProblemSetter of type $(typeof(pset))."
    pset
end
