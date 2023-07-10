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

function ProblemUpdater(par_getter::AbstractProblemParGetter, u0_keys, p_keys)
    ProblemUpdater(
        par_getter, ProblemParSetter(Axis(u0_keys),Axis(p_keys),Axis(keys(par_getter))))
end

par_setter(pu::ProblemUpdater) = pu.pset
par_getter(pu::ProblemUpdater) = pu.pget

function (pu::ProblemUpdater)(prob)
    x = pu.pget(pu, prob)
    update_statepar(pu.pset, x, prob)
end

"AbstractProblemUpdater that returns the original ODEProblem."
struct NullProblemUpdater{PG <: AbstractProblemParGetter, PS <: AbstractProblemParSetter} <: AbstractProblemUpdater; end
(pu::NullProblemUpdater)(prob) = prob


"""
    KeysProblemParGetter(source_keys::NTuple{N,Symbol})

Provices callable `(pg::KeysProblemParGetter)(pu::ProblemUpdater, prob)]`.    
To be used to get the parameters/state vector to be set by `ProblemUpdater`.

Initialize with an NTuple of symbols that index into 
`vcat(label_state(pu.pset, prob.u0), label_par(pu.pset, prob.p))`.
"""
struct KeysProblemParGetter{N} <: AbstractProblemParGetter
    source_keys::NTuple{N,Symbol}
    dest_keys::NTuple{N,Symbol}
    # type parameter already enfources same length
    # KeysProblemParGetter(source::NTuple{N,Symbol},dest::NTuple{N,Symbol}) 
    #     length(source) == length(dest) || error(
    #         "KeysProblemParGetter expected same length of source and dest keys," *
    #         "but was supplied with source=("*join(source,",")*") and dest=("*join(dest,",")*")")
    #     new{N}(source,dest)
    # end
end
# also implement the keys methods of a new AbstractProblemParGetter
# so that it can be used by constructing a ProblemUpdater
Base.keys(pg::KeysProblemParGetter) = pg.dest_keys

function (pg::KeysProblemParGetter)(pu::ProblemUpdater, prob) 
    p = vcat(label_state(pu.pset, prob.u0), label_par(pu.pset, prob.p))
    SVector(getproperty.(Ref(p), pg.source_keys))
end

"Custom AbstractProblemParGetter used in tests that computes parameters to set."
struct DummyParGetter <: AbstractProblemParGetter; end
function (pg::DummyParGetter)(pu::ProblemUpdater, prob) 
    p = label_par(pu.pset, prob.p)
    #p_source = getproperty.(Ref(p), SA[:k_L])
    #(:k_R,:k_P)
    vcat(p.k_L, p.k_L[1]*10)
end    
Base.keys(pg::DummyParGetter) = (:k_R,:k_P)