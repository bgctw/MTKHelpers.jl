abstract type AbstractProblemUpdater end
abstract type AbstractProblemParGetter end

"""
    ProblemUpdater(AbstractProblemParGetter, AbstractProblemParSetter) 

Encapsulates updating an ODEProblem based on the problem itself by 
Callable `(pu::ProblemUpdater)(prob)`.

Must be initialized with a callable `AbstractProblemParGetter`, 
e.g. [`KeysProblemParGetter`](@ref)
and on a `AbstractProblemParSetter`,
e.g. [`ProblemParSetter`](@ref).
"""
struct ProblemUpdater{PG <: AbstractProblemParGetter, PS <: AbstractProblemParSetter} <: AbstractProblemUpdater
    pget::PG
    pset::PS
end

function (pu::ProblemUpdater)(prob)
    x = pu.pget(pu, prob)
    update_statepar(pu.pset, x, prob)
end


"""
    KeysProblemParGetter(source_keys::NTuple{N,Symbol})

Provices callable `(pg::KeysProblemParGetter)(pu::ProblemUpdater, prob)]`.    
To be used to get the parameters/state vector to be set by `ProblemUpdater`.

Initialize with an NTuple of symbols that index into 
`vcat(label_state(pu.pset, prob.u0), label_par(pu.pset, prob.p))`.
"""
struct KeysProblemParGetter{N} <: AbstractProblemParGetter
    source_keys::NTuple{N,Symbol}
end

function (pg::KeysProblemParGetter)(pu::ProblemUpdater, prob) 
    p = vcat(label_state(pu.pset, prob.u0), label_par(pu.pset, prob.p))
    SVector(getproperty.(Ref(p), pg.source_keys))
end

"Custom AbstractProblemParGetter used in tests that computes parameters to set."
struct DummyParGetter <: AbstractProblemParGetter; end
function (pg::DummyParGetter)(pu::ProblemUpdater, prob) 
    p = label_par(pu.pset, prob.p)
    p_source = getproperty.(Ref(p), SA[:k_L])
    vcat(p_source, p_source[1]*10)
end    
