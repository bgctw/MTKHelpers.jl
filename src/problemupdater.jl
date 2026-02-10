abstract type AbstractProblemUpdater end

function get_concrete(pu::AbstractProblemUpdater)
    !isconcrete(pu) &&
        @warn "no concrete type implemented for ProblemUpdater of type $(typeof(pu))."
    pu
end

"""
Supertype for callables that implement 
    `(::AbstractProblemParGetter)(problem) -> updated_problem`

Concrete subtypes should implement function `keys`, so that an appropriate 
`AbstractODEProblemParSetter` can be constructed for [`ProblemUpdater`](@ref).
"""
abstract type AbstractProblemParGetter end

"""
    ProblemUpdater(par_getter, par_setter) 

Encapsulates updating an AbstractODEProblem based on the problem itself by 
Callable `(pu::ProblemUpdater)(prob)`.

Must be initialized with a callable `AbstractProblemParGetter`, 
e.g. [`KeysProblemParGetter`](@ref)
and on a `AbstractODEProblemParSetter`,
e.g. [`ODEProblemParSetter`](@ref).

There are special functions to construct ProblemUpdater based on
a given Problem:
- [`get_ode_problemupdater`](@ref)
"""
struct ProblemUpdater <: AbstractProblemUpdater
    pget::AbstractProblemParGetter
    pset::AbstractProblemParSetter
end

struct ProblemUpdaterConcrete{
    PG <: AbstractProblemParGetter,
    PS <: AbstractODEProblemParSetter,
} <: AbstractProblemUpdater
    pget::PG
    pset::PS
end

function get_concrete(pu::ProblemUpdater)
    ProblemUpdaterConcrete(par_getter(pu), get_concrete(par_setter(pu)))
end

ProblemUpdaterU = Union{ProblemUpdater, ProblemUpdaterConcrete}

"""
    get_ode_problemupdater(par_getter::AbstractProblemParGetter, u0, p)
    get_ode_problemupdater(par_getter::AbstractProblemParGetter, sys::AbstractSystem)

Construct a `ProblemUpdater` based on an constructed `ODEProblemParSetterConcrete`.     
"""
function get_ode_problemupdater(par_getter::AbstractProblemParGetter,
        sys::AbstractSystem)
    ProblemUpdater(par_getter, ODEProblemParSetter(sys, keys(par_getter)))
end,
function get_ode_problemupdater(par_getter::AbstractProblemParGetter, u0, p)
    error("deprecated: need a system to construct a ProblemPatSetter.")
    #ProblemUpdater(par_getter, ODEProblemParSetter(u0, p, keys(par_getter)))
end

@deprecate ProblemUpdater(par_getter, u0_keys, p_keys) get_ode_problemupdater(par_getter,
    u0_keys, p_keys)

par_setter(pu::ProblemUpdaterU) = pu.pset
par_getter(pu::ProblemUpdaterU) = pu.pget

function (pu::ProblemUpdaterU)(prob)
    x = pu.pget(pu, prob)
    remake(prob, x, pu.pset)
end

"""
AbstractProblemUpdater that returns the original AbstractODEProblem.
"""
struct NullProblemUpdater <: AbstractProblemUpdater end
(pu::NullProblemUpdater)(prob) = prob

isconcrete(::NullProblemUpdater) = true

"""
    KeysProblemParGetter(mapping::NTuple{N,Pair{Symbol, Symbol}, keys_state)

Provides callable `(pg::KeysProblemParGetter)(pu::ProblemUpdater, prob), keys_state]`.    
To be used to get the parameters/state vector to be set by `ProblemUpdater`.

Initialize with an mapping of NTuples of symbols (source -> target) that index into 
either `get_state_labeled(pu.pset, prob)` or `get_par_labeled(pu.pset, prob))`.
Argument `keys_state` is a Tuple or Vector that iterates the Symbols in the state of an 
ODEProblem. It is required to know from which part of the problem to extract.
"""
struct KeysProblemParGetter{N} <: AbstractProblemParGetter
    source_keys::NTuple{N, Symbol}
    dest_keys::NTuple{N, Symbol}
    is_in_state::SVector{N, Bool}
    # type parameter already enfources same length
end

function KeysProblemParGetter(mapping::NTuple{N, Pair{Symbol, Symbol}}, prob::AbstractODEProblem) where {N}
    sys = get_system(prob)
    KeysProblemParGetter(mapping, sys)
end

function KeysProblemParGetter(mapping::NTuple{N, Pair{Symbol, Symbol}}, sys::AbstractSystem) where {N}
    keys_state = symbol_op.(unknowns(sys))
    KeysProblemParGetter(mapping, keys_state)
end

function KeysProblemParGetter(mapping::NTuple{N, Pair{Symbol, Symbol}},
        keys_state::Union{AbstractVector{Symbol}, NTuple{NU, Symbol}}) where {N, NU}
    source_keys = ntuple(i -> first(mapping[i]), N)
    dest_keys = ntuple(i -> last(mapping[i]), N)
    if length(unique(dest_keys)) != length(dest_keys)
        error("Expected destination of pattern to be distinct, but was not.")
    end
    is_in_state = SVector{N}(k ∈ keys_state for k in source_keys)
    KeysProblemParGetter(source_keys, dest_keys, is_in_state)
end

# also implement the keys methods of a new AbstractProblemParGetter
# so that it can be used by constructing a ProblemUpdater
Base.keys(pg::KeysProblemParGetter) = pg.dest_keys

function (pg::KeysProblemParGetter{N})(pu::AbstractProblemUpdater, prob) where {N}
    pset = pu.pset
    u0l = get_state_labeled(pset, prob)
    pl = get_par_labeled(pu.pset, prob)
    T = promote_type(eltype(u0l), eltype(pl))
    #u0p = vcat(u0l, pl) # inferred Any with runtime dispatch []
    #vcat((u0p[k] for k in pg.source_keys)...)::Vector{T}
    #ftmp = (k) -> k ∈ keys(u0l) ? u0l[k] : pl[k]
    #vcat((ftmp(k) for k in pg.source_keys)...)::Vector{T}
    #Main.@infiltrate_main
    ftmp = (i, k) -> pg.is_in_state[i] ? u0l[k] : pl[k]
    vcat((ftmp(i, k) for (i, k) in enumerate(pg.source_keys))...)::Vector{T}
end

"""
Custom AbstractProblemParGetter used in tests that computes parameters to set.
update k_R = k_L and m = k_L*10
"""
struct DummyParGetter <: AbstractProblemParGetter end
function (pg::DummyParGetter)(pu::AbstractProblemUpdater, prob)
    p = get_par_labeled(pu.pset, prob)
    vcat(p[:k_L], p[:k_L] * 10)
end
Base.keys(pg::DummyParGetter) = (:k_R, :m)
