function is_symbolicarray(s::SymbolicUtils.BasicSymbolic)
    #istree(s) && (operation(s) == getindex) 
    iscall(s) && (operation(s) == getindex)
end

"""
    pos_of_base_nums(st)

Collect all occurrences of base_nums in a sequence of BasicSymbolics
`st = vcat(unknowns(sys), parameters(sys))`.
Return a `Dict{SymbolicUtils.BasicSymbolic,Vector{Int}}`.
"""
function pos_of_base_nums(st)
    dpos = Dict{SymbolicUtils.BasicSymbolic,Vector{Int}}()
    #(pos, sti) = first(enumerate(st))
    for (pos, sti) in enumerate(st)
        _bnum = base_num(sti)
        _posvec = get(dpos, _bnum, Vector{Int}())
        if isempty(_posvec)
            dpos[_bnum] = _posvec
        end
        push!(_posvec, pos)
    end
    dpos
end

"""
    base_num(s)

Get the symbol without an index, e.g. p[1] -> p, or x(t)[1] -> x(t).   
"""
function base_num(s::SymbolicUtils.BasicSymbolic)
    #!istree(s) ? s :
    !iscall(s) ? s :
    #operation(s) == getindex ? base_num(first(arguments(s))) : base_num(operation(s)) 
    operation(s) == getindex ?
    base_num(first(arguments(s))) : s
end
base_num(s) = s

"""
    get_system_symbol_dict(sys::AbstractODESystem)

Construct a `Dict{Symbol => Num}` for all properties in `sys`.
"""
function get_system_symbol_dict(sys::AbstractODESystem)
    # if there are no observed, return type is Dict(Any,Any) -> need conditional
    dicts = (
        get_base_num_dict(unknowns(sys)),
        get_base_num_dict(parameters(sys)),
        get_scalarized_num_dict(unknowns(sys)),
        get_scalarized_num_dict(parameters(sys)),
        get_base_num_dict(getproperty.(observed(sys), :lhs)),
        get_scalarized_num_dict(getproperty.(observed(sys), :lhs))
    )
    # only merge nonemtpy dictionaries, otherwise the eltype becomes Any
    dicts_nonempty = filter(d -> !isempty(d), dicts)
    merge(dicts_nonempty...)
end

function get_system_symbol_dict(sys, cv;
    system_symbol_dict=get_system_symbol_dict(sys))
    ssd = system_symbol_dict
    pd = vcat((length(cv[k]) == 1 ?
               ssd[k] => cv[k] :
               Symbolics.scalarize(ssd[k] .=> collect(cv[k])) for k in keys(cv))...)
end

"""
    get_base_num_dict(nums)

Return a Dictionary of Symbol -> Num, for each unique `base_num.(nums)`
"""
function get_base_num_dict(nums, f_symbol=symbol_op)
    @chain nums begin
        base_num.()
        unique()
        f_symbol.(_) .=> _
        Dict()
    end
end

"""
Return a Dictionar y of Symbol -> Num but only for each symbolicarray among nums
"""
function get_scalarized_num_dict(nums)
    #nums_tree = filter(x -> istree(x) && (operation(x) == getindex), nums)
    nums_tree = filter(x -> iscall(x) && (operation(x) == getindex), nums)
    # need to take care of type, because empty case returns Dict{Any,Any}
    isempty(nums_tree) ? Dict{Symbol,eltype(nums)}() :
    # in Julia 1.6 nums maybe Any, but Dict has values of more specific type -> convert
    convert(Dict{Symbol,eltype(nums)},
        Dict(Symbol.(nums_tree) .=> nums_tree))::Dict{Symbol,eltype(nums)}
    # chain version is not inferred
    # @chain nums begin 
    #     filter(istree, _)
    #     filter(x -> , _)
    #     Symbol.(_) .=> _
    #     # Dict()        
    # end
end



# function get_system_symbol_dict(sys::AbstractSystem,
#     string_sys::String = string(nameof(sys)))
#     prefix = isempty(string_sys) ? "" : string_sys * "₊"
#     Dict(Symbol(prefix * string(p)) => getproperty(sys, p) for
#          p in propertynames(sys))
# end
# function get_system_symbol_dict(systems...)
#     dicts = map(systems) do sys
#         get_system_symbol_dict(sys)
#     end
#     merge(dicts...)
# end

@deprecate strip_deriv_num(x) symbol_op(x)

"""
    remake_cv(prob::AbstractODEProblem, paropt::ComponentVector; 
        num_dict_state = get_base_num_dict(unknowns(get_system(prob))),
        num_dict_par = get_base_num_dict(parameters(get_system(prob)))
    
Creates a new problem with components in `u0` and `p` begin updated, for problems
with an associated ODESystem.
For doing this more efficiently when repeating, pre-compute the `num_dict_state` and `num_dict_par` once,
and provide it to the function.

For discretized pde systems, some of the scalarized states are computed
by boundary conditions and are not in the state vector.
For those, supply argument `state_pos` giving indices of the states,
as found by [`get_1d_state_pos`](@ref).
"""
# function remake_cv(prob::AbstractODEProblem, paropt::ComponentVector;
#     state_pos=nothing,
#     num_dict_state=get_base_num_dict(unknowns(get_system(prob))),
#     num_dict_par=get_base_num_dict(parameters(get_system(prob))))
#     u0 = componentvector_to_numdict(paropt.state, num_dict_state; indices=state_pos)
#     p = componentvector_to_numdict(paropt.par, num_dict_par)
#     SciMLBase.remake(prob, u0=u0, p=p)
# end

"""
    system_num_dict(d, sys::AbstractSystem)
    system_num_dict(d, symbol_dict::AbstractDict)

Create a Dictionary Num=>value from symbolic Dictionary or ComponentVector.

Omit pairs where no Num was found.
"""
function system_num_dict(d, sys::AbstractSystem)
    system_num_dict(d, get_system_symbol_dict(sys))
end
# Dictionary or Vector of pairs Symbol -> number
function system_num_dict(d, symbol_dict::AbstractDict)
    Dict([get(symbol_dict, k, missing) => v
          for
          (k, v) in d if !ismissing(get(symbol_dict, k, missing))])
end
function system_num_dict(ca::CA.ComponentVector, symbol_dict::AbstractDict)
    componentvector_to_numdict(ca, symbol_dict)
end

function componentvector_to_numdict(cv::ComponentVector{T}, num_dict::Dict{Symbol,S};
    indices=nothing) where {T,S}
    #kcv = first(keys(cv)) # kcv=last(keys(cv))
    get_pairs = (kcv) -> begin # returns either Pair or Vector{Pair}
        num = num_dict[kcv]
        num_s = Symbolics.scalarize(num)
        # is_vector = is_symbolicarray(num) 
        # _pairs = !is_vector ? num => cv[symbol_op_scalar(num)] :
        #     isnothing(indices) ? 
        #     num_s .=> cv[symbol_op_scalar(num)] : 
        #     num_s[indices] .=> cv[symbol_op_scalar(num)]
        _pairs = isnothing(indices) ?
                 num_s .=> cv[symbol_op_scalar(num)] :
                 num_s[indices] .=> cv[symbol_op_scalar(num)]
    end
    _pairs_gen = (get_pairs(kcv) for kcv in keys(cv))
    _pairs_all = reduce(vcat, _pairs_gen; init=Vector{T}())
    isempty(_pairs_all) && return Dict{SymbolicUtils.BasicSymbolic{Real},eltype(cv)}()
    Dict(_pairs_all...)::Dict{SymbolicUtils.BasicSymbolic{Real}}
end
# fallback for empty subvectors of a ComponentArray
function componentvector_to_numdict(cv::SubArray{T}, num_dict::Dict{Symbol,S}) where {T,S}
    Dict{SymbolicUtils.BasicSymbolic{Real},T}()
end



"""
    expand_base_num(num, sys::AbstractSystem) 

Return the scalaized symbols for a Num of a PDESystem discretized along one dimension.
The result is subsetted by argument `state_pos`, which defaults to `get_1d_state_pos(sys)`.
"""
function expand_base_num(num, sys::AbstractSystem)
    state_pos = get_1d_state_pos(sys)
    expand_base_num(num)[state_pos]
end

# function expand_base_num(num, state_pos::AbstractVector{Int})
#     Symbolics.scalarize(num)[state_pos]
# end

"""
    symbol_num_getindex(num::SymbolicUtils.Symbolic)

return Symbol(num) with the exception of getindex, then instead of forms like
`Symbol("getindex(a(t), 1)")` return `Symbol("(a(t))[1]")`.
"""
function symbol_num_getindex(num::SymbolicUtils.Symbolic)
    iscall(num) && (operation(num) == getindex) ? Symbol(string(num)) : Symbol(num)
end

"""
    expand_base_num_axes(cv::ComponentVector, sys::AbstractSystem)

Change the axis of a componentvector to replace vector-valued entries by
their respective scalarized symbols. 
For array-symbols, assume only one index and order symbols by that index.
E.g. `Y=[1:14]` gets transformed to `[Y(t))[2], (Y(t))[3], ..., (Y(t)[15]]`,
where the `Y(t))[i]` are part of `unknows(sys)` and can appear in any order in
the system.
"""
function expand_base_num_axes(cv::ComponentVector, sys::AbstractSystem)
    @warn("deprecated: use expand_base_num_axes(cv, get_scalar_num_map(sys))")
    scalar_num_map = get_scalar_num_map(sys)
    expand_base_num_axes(cv, scalar_num_map)
end
expand_base_num_axes(cv::UnitRange, sys::AbstractSystem) = cv

function expand_base_num_axes(
    cv::ComponentVector, scalar_num_map::Dict{<:SymbolicUtils.BasicSymbolic})
    scalar_num_map_sym = Dict(symbol_op(k) => v for (k, v) in scalar_num_map)
    expand_base_num_axes(cv, scalar_num_map_sym)
end
function expand_base_num_axes(cv::ComponentVector, scalar_num_map::Dict{Symbol})
    #state_pos = get_1d_state_pos(sys)
    #k = first(keys(cv))
    #k = last(keys(cv))
    tmp = map(keys(cv)) do k
        length_x = length(cv[CA.KeepIndex(k)])
        haskey(scalar_num_map, k) ||
            error("Expected all keys in cv to map to numerics." *
                  "but found no match for '$k'.")
        syms_k = scalar_num_map[k]
        length_x == 1 && first(syms_k)
        length(syms_k) == length_x ||
            error("Expected component entry $k to hold $(length(syms_k)) entries, " *
                  "but was $(length_x)")
        symbol_num_getindex.(syms_k)
    end
    syms = reduce(vcat, tmp) # vcat of SA of Symbols is better than vcat of CA
    # SA transferred to type - not typestable but saves precompilation compared to
    # vcat of Sub-ComponentArrays
    ax_scalar = Axis(syms)
    attach_axis(getdata(cv), ax_scalar)
end
expand_base_num_axes(cv::UnitRange, scalar_num_map::Dict) = cv

"""
    get_scalar_num_map(sys::AbstractSystem)

Return a mapping of each uniuqe base_num of the system
to original scalarized BasicSymbolics Nums used in the System.
The order of the nums of symbolic arrays is ascending and differs 
from the order in the system.
The translation in the problem is then taken care of by `remake(Problem, popt, pset)`,
which uses this map to create the Symbolic mapping or indexing
into the state.

In order to get a mapping from symbols use [`get_scalar_num_map_sym`](@ref)
"""
function get_scalar_num_map(sys::AbstractSystem)
    st = vcat(unknowns(sys), parameters(sys))
    dpos = pos_of_base_nums(st)
    #dpos_sym = Dict(symbol_op(p.first) => p.second for p in dpos)
    dpos_sym = Dict(p.first => p.second for p in dpos)
    (; st, dpos_sym)
    scalar_nums_map = Dict(k => getindex.(Ref(st), dpos_sym[k]) for k in keys(dpos_sym))
end,
function get_scalar_num_map_sym(sys::AbstractSystem)
    scalar_num_map = get_scalar_num_map(sys)
    Dict(symbol_op(k) => v for (k, v) in scalar_num_map)
end


"""
    axis_of_nums(nums)

Return an axis with keys corresponding to each symbolic array and values
to the indices in the given vector of nums.
Assumes that symbolic arrays are consecutive positions in vector.
"""
function axis_of_nums(nums::AbstractVector)
    pos_nums = indices_of_nums(nums)
    pos_nums_scalar = [k => length(v) == 1 ? first(v) : v for (k, v) in pos_nums]
    #CA._component_axis(first(axes(CA.ComponentArray(;pos_nums_scalar...))))
    first(getaxes(CA.ComponentArray(; pos_nums_scalar...)))
end

function indices_of_nums(nums)
    op_syms = symbol_op.(nums)
    #vec_of_pairs = Vector{Pair{Symbol, Union{Int,UnitRange{Int}}}}()
    vec_of_pairs = Vector{Pair{Symbol,UnitRange{Int}}}()
    for (pos, sym) in enumerate(op_syms)
        if length(vec_of_pairs) == 0
            push!(vec_of_pairs, sym => pos:pos)
        else
            current_pair = vec_of_pairs[end]
            if sym != first(current_pair)
                push!(vec_of_pairs, sym => pos:pos)
            else
                current_range = Base.last(current_pair)
                vec_of_pairs[end] = sym => first(current_range):pos
            end
        end
    end
    return vec_of_pairs
end


