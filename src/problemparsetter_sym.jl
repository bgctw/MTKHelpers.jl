struct ProblemParSetter_sym{NS,NP,NO,IT} <: AbstractProblemParSetter
    optinfo::NTuple{NO, Tuple{Symbol,Symbol,IT}}
    statemap::NTuple{NS, IT} # position_in_u0 -> posopt, 0: not optimized
    parmap::NTuple{NP, IT} # position_in_p -> posopt, 0: not optimized
    statesyms::NTuple{NS, Symbol}
    parsyms::NTuple{NP, Symbol}
    # isuopt::NTuple{NP, Bool}
    # ispopt::NTuple{NP, Bool}
end

"""
    ProblemParSetter_sym(state_names,par_names,popt_names) 
    ProblemParSetter_sym(sys::ODESystem, popt_names; strip=false) 

Helps keeping track of a subset of initial states and paraemters to be optimized.

# Arguments
- `state_names`: Tuple or AbstractVector of all the initial states of the problem
- `par_names`: all the parameters of the problem
- `popt_names`: the parameter/initial states to be optimized.

The states and parameters can be extracted from an `ModelingToolkit.ODESystem`.
If `strip=true`, then namespaces of parameters of a composed system are removed, 
e.g. `subcomp₊p` becomes `p`.
"""
function ProblemParSetter_sym(
    state_syms::NTuple{NS,Symbol},par_syms::NTuple{NP,Symbol},popt_syms::NTuple{NO,Symbol}, 
    ::Val{allow_missing_popt}=Val(true); it::Type{IT}=Int64) where {NO, NS,NP,IT,allow_missing_popt}
    posu_orig = Dict(parname => convert(IT,pos) for (pos, parname) in enumerate(state_syms))
    posp_orig = Dict(parname => convert(IT,pos) for (pos, parname) in enumerate(par_syms))
    #
    optinfo_u = Tuple((:state, statesym, posu_orig[statesym]) for 
        statesym in popt_syms if haskey(posu_orig, statesym))
    optinfo_p = Tuple((:par, parsym, posp_orig[parsym]) for 
        parsym in popt_syms if haskey(posp_orig, parsym))
    optinfo = (optinfo_u..., optinfo_p...)
    #    
    map_posu_to_posopt = Dict(pos => posopt for 
        (posopt, (source, parname, pos)) in enumerate(optinfo) if source == :state)
    map_posp_to_posopt = Dict(pos => posopt for 
        (posopt, (source, parname, pos)) in enumerate(optinfo) if source == :par)
        # 
    statemap = NTuple{NS}((get(map_posu_to_posopt,posu,0) for posu in one(IT):NS))
    parmap = NTuple{NP}((get(map_posp_to_posopt,posp,0) for posp in one(IT):NP))
    # isuopt = NTuple{NS}((i != 0 for i in statemap))
    # ispopt = NTuple{NP}((i != 0 for i in parmap))
    # @show NS, NP, NO
    # @show optinfo
    # @show statemap
    # @show parmap
    popt_syms_in = popt_syms
    if allow_missing_popt && (NO != length(optinfo))
        optinfo_psyms = getindex.(optinfo,2)
        miss = setdiff(popt_syms, optinfo_psyms)
        @warn(
            "missing optimization parameters in system: " * join(miss,", "))
        popt_syms_in = filter((x -> !isnothing(findfirst(==(x), optinfo_psyms))), popt_syms)
        ps = ProblemParSetter_sym{NS,NP,length(optinfo),IT}(
            optinfo, statemap, parmap, state_syms, par_syms,)     
    else
        # note using NO here for type stability
        ps = ProblemParSetter_sym{NS,NP,NO,IT}(
            optinfo, statemap, parmap, state_syms, par_syms,)     
    end
    # make sure that states go first: smallest p position > highest u position            
    symbols_paropt(ps) == popt_syms_in || error(
        "Expected states to optimize before parameters to optimize. But got $popt_syms_in.")
    ps
end

function ProblemParSetter_sym(state_names,par_names,popt_names) 
    # not type-stable
    ProblemParSetter_sym(
        Tuple(i for i in symbol.(state_names)),
        Tuple(i for i in symbol.(par_names)),
        Tuple(i for i in symbol.(popt_names))
    )
end

function ProblemParSetter_sym(sys::ODESystem,popt_names; strip=false) 
    ft = strip ? strip_namespace : identity
    state_names = ft.(symbol.(states(sys)))
    par_names = ft.(symbol.(parameters(sys)))
    ProblemParSetter_sym(state_names, par_names, popt_names)
end

function count_state(::ProblemParSetter_sym{NS}) where {NS}; NS; end,
function count_par(::ProblemParSetter_sym{NS,NP}) where {NS,NP}; NP; end,
function count_paropt(::ProblemParSetter_sym{NS,NP,NO}) where {NS,NP,NO}; NO; end

@deprecate count_states(ps::ProblemParSetter_sym)  count_state(ps)


@deprecate statesyms(ps::ProblemParSetter_sym) symbols_state(ps)
@deprecate parsyms(ps::ProblemParSetter_sym) symbols_par(ps)
@deprecate paroptsyms(ps::ProblemParSetter_sym) symbols_paropt(ps)

function symbols_state(ps::ProblemParSetter_sym); ps.statesyms; end,
function symbols_par(ps::ProblemParSetter_sym); ps.parsyms; end,
function symbols_paropt(ps::ProblemParSetter_sym)
    first.(Base.Iterators.drop.(ps.optinfo,1))
end

function update_statepar(ps::ProblemParSetter_sym, popt, prob::ODEProblem) 
    u0,p = update_statepar(ps, popt, prob.u0, prob.p)
    remake(prob; u0, p)
end

function update_statepar(ps::ProblemParSetter_sym, popt, u0::TU, p::TP) where {TU,TP}
    @assert length(popt) == count_paropt(ps) "expected $(count_paropt(ps)) parameters for "
    # @show TU
    # u2 = TU((ps.statemap[i] == 0 ? u0[i] : popt[ps.statemap[i]] for i in 1:length(u0)))
    # p2 = TP((ps.parmap[i] == 0 ? p[i] : popt[ps.parmap[i]] for i in 1:length(p)))
    # (u2,p2)
    # take care that popt might be of different type, like FowwardDiff
    # need to convert to new type
    u00 = map(x -> x * zero(eltype(popt)), u0)
    #eachindex(u0) does not work it iterates keys instead of integer indices
    u0g = (ps.statemap[i] == 0 ? u0[i] : popt[ps.statemap[i]] for i in 1:length(u0))
    #u0new = typed_from_generator(TU, u00, u0g)
    u0new = typed_from_generator(u00, u0g)
    #u0new = typeof(u00)(u0g)
    #
    p0 = map(x -> x * zero(eltype(popt)), p)
    # pg = (ps.parmap[i] == 0 ? p[i] : popt[ps.parmap[i]] for i in 1:length(p))
    # pnew = convert(typeof(p0), p0 .+ pg)::typeof(p0)
    pg = (ps.parmap[i] == 0 ? p[i] : popt[ps.parmap[i]] for i in 1:length(p))
    #pnew = typed_from_generator(TP, p0, pg)
    pnew = typed_from_generator(p0, pg)
    #pnew = typeof(p0)(pg)
    #
    (u0new, pnew)
end

# function typed_from_generator(type::Type, v0, vgen) where T 
#     if type <: AbstractArray
#     end
#     typeof(v0)(vgen)
# end
# # AbstractVector{T} is not more specific than ::Type, need to support all used concrete
# #typed_from_generator(::Type{AbstractVector{T}}, v0, vgen) where T = convert(typeof(v0), v0 .+ vgen)::typeof(v0)
# typed_from_generator(::Type{Vector{T}}, v0, vgen) where T = convert(typeof(v0), v0 .+ vgen)::typeof(v0)


typed_from_generator(v0, vgen) = typeof(v0)(vgen)
function typed_from_generator(v0::AbstractVector, vgen) 
    # convert to std vector, because typeof(v0) does not contain all entries
    T = Vector{eltype(v0)}  
    convert(T, v0 .+ vgen)::T
end

function get_paropt(ps::ProblemParSetter_sym, prob::ODEProblem; kwargs...)
    get_paropt(ps, prob.u0, prob.p; kwargs...)
end
function get_paropt_labeled(ps::ProblemParSetter_sym, prob::ODEProblem; kwargs...)
    get_paropt_labeled(ps, prob.u0, prob.p; kwargs...)
end

function get_paropt(ps::ProblemParSetter_sym, u0::AbstractVector, p::AbstractVector) 
    v = [(first(t) == :par) ? p[last(t)] : u0[last(t)] for t in ps.optinfo]
end

function get_paropt_labeled(ps::ProblemParSetter_sym, u0::AbstractVector, p::AbstractVector) 
    v = get_paropt(ps, u0, p)
    label_paropt(ps,v)
end


function get_paropt(ps::ProblemParSetter_sym{NS,NP,NO}, u0, p) where {NS,NP,NO}
    t0 = NTuple{NO}(((first(t) == :par) ? p[last(t)] : u0[last(t)] 
        for t in ps.optinfo))
    # need to explicitly assure full type for julia 1.6 for type stability
    t1 = t0::NTuple{NO,eltype(t0)} 
end

function get_paropt_labeled(ps::ProblemParSetter_sym{NS,NP,NO}, u0, p) where {NS,NP,NO}
    t1 = get_paropt(ps,u0,p)
    label_paropt(ps,t1)
end

"""
    label_state(ps::ProblemParSetter_sym, u::AbstractVector)
    label_par(ps::ProblemParSetter_sym, par::AbstractVector)
    label_paropt(ps::ProblemParSetter_sym, popt::AbstractVector)

Produce a labeled version of a sequence of initial states, parameters, or
optimized parameters respectively.
The return type differs given the input
- SVector -> SLVector
- NTuple -> NamedTuple
- AbstractVector -> LArray
"""
function label_state(ps::ProblemParSetter_sym, u::AbstractVector); LArray{symbols_state(ps)}(u); end,
function label_par(ps::ProblemParSetter_sym, par::AbstractVector); LArray{symbols_par(ps)}(par); end,
function label_paropt(ps::ProblemParSetter_sym, popt::AbstractVector); LArray{symbols_paropt(ps)}(popt); end

label_state(ps::ProblemParSetter_sym, u::SVector) = SLVector(label_state(ps, Tuple(u)))
label_state(ps::ProblemParSetter_sym, u::NTuple) = NamedTuple{symbols_state(ps)}(u)

label_par(ps::ProblemParSetter_sym, par::SVector) = SLVector(label_par(ps, Tuple(par)))
label_par(ps::ProblemParSetter_sym, par::NTuple) = NamedTuple{symbols_par(ps)}(par)

label_paropt(ps::ProblemParSetter_sym, popt::SVector) = SLVector(label_paropt(ps, Tuple(popt)))
label_paropt(ps::ProblemParSetter_sym, popt::NTuple) = NamedTuple{symbols_paropt(ps)}(popt)

# extends Base.merge to work on SVector
merge(x::T, y::NamedTuple) where T<:SLArray = T(merge(NamedTuple(x),y)...)

# and on Labelled Arrays
function merge(x::T, y::NamedTuple) where T<:LArray
    xnew = deepcopy(x)
    for (key,val) in pairs(y)
        xnew[key] = val
    end
    xnew
end

name_state(ps::ProblemParSetter_sym, state::AbstractVector) = NamedArray(state, (collect(symbols_state(ps)),))
name_par(ps::ProblemParSetter_sym, par::AbstractVector) = NamedArray(par, (collect(symbols_par(ps)),))
name_paropt(ps::ProblemParSetter_sym, paropt::AbstractVector) = NamedArray(paropt, (collect(symbols_paropt(ps)),))



