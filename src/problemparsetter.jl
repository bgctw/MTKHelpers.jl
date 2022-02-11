struct ProblemParSetter{NS,NP,NO,IT}
    optinfo::NTuple{NO, Tuple{Symbol,Symbol,IT}}
    statemap::NTuple{NS, IT} # position_in_u0 -> posopt, 0: not optimized
    parmap::NTuple{NP, IT} # position_in_p -> posopt, 0: not optimized
    statesyms::NTuple{NS, Symbol}
    parsyms::NTuple{NP, Symbol}
    # isuopt::NTuple{NP, Bool}
    # ispopt::NTuple{NP, Bool}
end

"""
    ProblemParSetter(state_names,par_names,popt_names) 
    ProblemParSetter(sys::ODESystem, popt_names; strip=false) 

Helps keeping track of a subset of initial states and paraemters to be optimized.

# Arguments
- `state_names`: Tuple or Vector of all the initial states of the problem
- `par_names`: all the parameters of the problem
- `popt_names`: the parameter/initial states to be optimized.

The states and parameters can be extracted from an `ModelingToolkit.ODESystem`.
If `strip=true`, then namespaces of parameteres of a composed system are removed, 
e.g. `subcomp₊p` becomes `p`.
"""
function ProblemParSetter(
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
    if allow_missing_popt
        NO != length(optinfo) && @warn(
        "missing optimization parameters in system: " * join(
        setdiff(par_syms, getindex.(optinfo,2))
        ,", "))
        ps = ProblemParSetter{NS,NP,length(optinfo),IT}(
            optinfo, statemap, parmap, state_syms, par_syms,)     
    else
        # note using NO here for type stability
        ps = ProblemParSetter{NS,NP,NO,IT}(
            optinfo, statemap, parmap, state_syms, par_syms,)     
    end
end

function ProblemParSetter(state_names,par_names,popt_names) 
    # not type-stable
    ProblemParSetter(
        Tuple(i for i in symbol.(state_names)),
        Tuple(i for i in symbol.(par_names)),
        Tuple(i for i in symbol.(popt_names))
    )
end

function ProblemParSetter(sys::ODESystem,popt_names; strip=false) 
    state_names = strip ? strip_namespace.(statesyms(sys)) : statesyms(sys)
    par_names = strip ? strip_namespace.(parsyms(sys)) : parsyms(sys)
    ProblemParSetter(state_names, par_names, popt_names)
end

"""
    count_states(::ProblemParSetter) 
    count_par(::ProblemParSetter) 
    count_paropt(::ProblemParSetter) 

Report the number of problem states, problem parameters and optimized parameters
respectively.    
"""
function count_states(::ProblemParSetter{NS}) where {NS}; NS; end,
function count_par(::ProblemParSetter{NS,NP}) where {NS,NP}; NP; end,
function count_paropt(::ProblemParSetter{NS,NP,NO}) where {NS,NP,NO}; NO; end


"""
    statesyms(ps::ProblemParSetter)
    parsyms(ps::ProblemParSetter)
    paroptsyms(ps::ProblemParSetter)

Report the names, i.e. symbols of problem states, problem parameters and 
optimized parameters respectively.    
"""
function statesyms(ps::ProblemParSetter); ps.statesyms; end,
function parsyms(ps::ProblemParSetter); ps.parsyms; end,
function paroptsyms(ps::ProblemParSetter)
    first.(Base.Iterators.drop.(ps.optinfo,1))
end

"""
    prob_new = update_statepar(ps::ProblemParSetter, popt, prob::ODEProblem) 
    u0new, pnew = update_statepar(ps::ProblemParSetter, popt, u0, p) 

Return an updated problem or updates states and parameters where
values corresponding to positions in `popt` are set.
"""
function update_statepar(ps::ProblemParSetter, popt, prob::ODEProblem) 
    u0,p = update_statepar(ps, popt, prob.u0, prob.p)
    remake(prob; u0, p)
end

function update_statepar(ps::ProblemParSetter, popt, u0::TU, p::TP) where {TU,TP}
    @assert length(popt) == count_paropt(ps) "expected $(count_paropt(ps)) parameters for "
    # @show TU
    # u2 = TU((ps.statemap[i] == 0 ? u0[i] : popt[ps.statemap[i]] for i in 1:length(u0)))
    # p2 = TP((ps.parmap[i] == 0 ? p[i] : popt[ps.parmap[i]] for i in 1:length(p)))
    # (u2,p2)
    # take care that popt might be of different type, like FowwardDiff
    # need to convert to new type
    u00 = map(x -> x * zero(eltype(popt)), u0)
    u0g = (ps.statemap[i] == 0 ? u0[i] : popt[ps.statemap[i]] for i in 1:length(u0))
    u0new = typeof(u00)(u0g)
    #
    p0 = map(x -> x * zero(eltype(popt)), p)
    # pg = (ps.parmap[i] == 0 ? p[i] : popt[ps.parmap[i]] for i in 1:length(p))
    # pnew = convert(typeof(p0), p0 .+ pg)::typeof(p0)
    pg = (ps.parmap[i] == 0 ? p[i] : popt[ps.parmap[i]] for i in 1:length(p))
    pnew = typeof(p0)(pg)
    #
    (u0new, pnew)
end

function update_statepar(ps::ProblemParSetter, popt, u0::Vector{T}, p::Vector{T}) where T
    @assert length(popt) == count_paropt(ps) "expected $(count_paropt(ps)) parameters for "
    #u2 = [ps.statemap[i] == 0 ? u0[i] : popt[ps.statemap[i]] for i in 1:length(u0)]
    #p2 = [ps.parmap[i] == 0 ? p[i] : popt[ps.parmap[i]] for i in 1:length(p)]
    # take care that popt might be of different type, like FowwardDiff
    # need to convert to new type
    u00 = map(x -> x * zero(eltype(popt)), u0)
    u0g = (ps.statemap[i] == 0 ? u0[i] : popt[ps.statemap[i]] for i in axes(u0,1))
    u0new = convert(typeof(u00), u00 .+ u0g)::typeof(u00)
    #
    p0 = map(x -> x * zero(eltype(popt)), p)
    pg = (ps.parmap[i] == 0 ? p[i] : popt[ps.parmap[i]] for i in axes(p,1))
    pnew = convert(typeof(p0), p0 .+ pg)::typeof(p0)
    #
    (u0new, pnew)
end


"""
    get_paropt(ps::ProblemParSetter, prob::ODEProblem; kwargs...)
    get_paropt(ps::ProblemParSetter, u0, p)

    get_paropt_labeled(ps::ProblemParSetter, prob::ODEProblem; kwargs...)
    get_paropt_labeled(ps::ProblemParSetter, u0, p)

Extract the initial states and parameters corresponding to the positions
that are optimized.    
If both u0 and p are vectors, the result is a vector, otherwise the result
is a Tuple.

The _lebeled versions additionally call `label_paropt` (see [`label_state`](@ref)) 
on the return value.
"""
function get_paropt(ps::ProblemParSetter, prob::ODEProblem; kwargs...)
    get_paropt(ps, prob.u0, prob.p; kwargs...)
end
function get_paropt_labeled(ps::ProblemParSetter, prob::ODEProblem; kwargs...)
    get_paropt_labeled(ps, prob.u0, prob.p; kwargs...)
end

function get_paropt(ps::ProblemParSetter, u0::Vector, p::Vector) 
    v = [(first(t) == :par) ? p[last(t)] : u0[last(t)] for t in ps.optinfo]
end

function get_paropt_labeled(ps::ProblemParSetter, u0::Vector, p::Vector) 
    v = get_paropt(ps, u0, p)
    label_paropt(ps,v)
end


function get_paropt(ps::ProblemParSetter{NS,NP,NO}, u0, p) where {NS,NP,NO}
    t0 = NTuple{NO}(((first(t) == :par) ? p[last(t)] : u0[last(t)] 
        for t in ps.optinfo))
    # need to explicitly assure full type for julia 1.6 for type stability
    t1 = t0::NTuple{NO,eltype(t0)} 
end

function get_paropt_labeled(ps::ProblemParSetter{NS,NP,NO}, u0, p) where {NS,NP,NO}
    t1 = get_paropt(ps,u0,p)
    label_paropt(ps,t1)
end

"""
    label_state(ps, u::Vector) = LArray{statesyms(ps)}(u)
    label_par(ps, par::Vector) = LArray{parsyms(ps)}(par)
    label_paropt(ps, popt::Vector) = LArray{paroptsyms(ps)}(popt)

Produce a labeled version of a sequence of initial states, parameters, or
optimized parameters respectively.
The return type differs given the input
- SVector -> SLVector
- NTuple -> NamedTuple
- Vector -> LArray
"""
function label_state(ps, u::Vector); LArray{statesyms(ps)}(u); end,
function label_par(ps, par::Vector); LArray{parsyms(ps)}(par); end,
function label_paropt(ps, popt::Vector); LArray{paroptsyms(ps)}(popt); end

label_state(ps, u::SVector) = SLVector(label_state(ps, Tuple(u)))
label_state(ps, u::NTuple) = NamedTuple{statesyms(ps)}(u)

label_par(ps, par::SVector) = SLVector(label_par(ps, Tuple(par)))
label_par(ps, par::NTuple) = NamedTuple{parsyms(ps)}(par)

label_paropt(ps, popt::SVector) = SLVector(label_paropt(ps, Tuple(popt)))
label_paropt(ps, popt::NTuple) = NamedTuple{paroptsyms(ps)}(popt)

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



