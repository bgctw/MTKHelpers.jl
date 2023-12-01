"""
    base_num(s)

Get the symbol without an index, e.g. p[1] -> p, or x(t)[1] -> x(t).   
"""
function base_num(s::SymbolicUtils.BasicSymbolic)
    !istree(s) ? s :
    #operation(s) == getindex ? base_num(first(arguments(s))) : base_num(operation(s)) 
    operation(s) == getindex ?
    base_num(first(arguments(s))) : s
end
base_num(s) = s

"""
    get_base_num_dict(nums)

Return a Dictionary of Symbol -> Num, for each unique `base_num.(nums)`
"""
function get_base_num_dict(nums)
    @chain nums begin
        base_num.()
        unique()
        symbol_op.(_) .=> _
        Dict()
    end
end

"""
    get_system_symbol_dict(sys::AbstractSystem, string_sys::String=string(nameof(sys)))
    get_system_symbol_dict(systems...)

Construct a `Dict{Symbol => Num}` for all properties in `sys`.
All Symbols are prefixed with `<string_sys>₊`

The second variant merges the dictionaries obtained from several systems.
"""
function get_system_symbol_dict(sys::AbstractODESystem)
    merge(
        # get_base_num_dict(ModelingToolkit.namespace_variables(sys)),
        # get_base_num_dict(ModelingToolkit.namespace_parameters(sys)),
        get_base_num_dict(states(sys)),
        get_base_num_dict(parameters(sys)))
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

function componentvector_to_numdict(cv::ComponentVector{T}, num_dict::Dict{Symbol, S};
        indices = nothing) where {T, S}
    tmp = @chain keys(cv) begin
        #filter(k -> k ∈ keys(num_dict), _)
        isnothing(indices) ?
        (Symbolics.scalarize(k) .=> cv[symbol_op(k)] for k in getindex.(Ref(num_dict), _)) :
        (Symbolics.scalarize(k)[indices] .=> cv[symbol_op(k)]
         for k in getindex.(Ref(num_dict), _))
        vcat(_...)
        Dict(_)
    end
    isempty(tmp) && return Dict{SymbolicUtils.BasicSymbolic{Real}, T}()
    tmp::Dict{SymbolicUtils.BasicSymbolic{Real}, T}
end
# fallback for empty subvectors of a ComponentArray
function componentvector_to_numdict(cv::SubArray{T}, num_dict::Dict{Symbol, S}) where {T, S}
    Dict{SymbolicUtils.BasicSymbolic{Real}, T}()
end

"""
    remake_cv(prob::AbstractODEProblem, paropt::ComponentVector; 
        num_dict_state = get_base_num_dict(states(get_system(prob))),
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
function remake_cv(prob::AbstractODEProblem, paropt::ComponentVector;
        state_pos = nothing,
        num_dict_state = get_base_num_dict(states(get_system(prob))),
        num_dict_par = get_base_num_dict(parameters(get_system(prob))))
    u0 = componentvector_to_numdict(paropt.state, num_dict_state; indices = state_pos)
    p = componentvector_to_numdict(paropt.par, num_dict_par)
    SciMLBase.remake(prob, u0 = u0, p = p)
end

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
    Dict([get(symbol_dict, k, missing) => v for
          (k, v) in d if !ismissing(get(symbol_dict, k, missing))])
end
function system_num_dict(ca::CA.ComponentVector, symbol_dict::AbstractDict)
    componentvector_to_numdict(ca, symbol_dict)
end

function indices_of_nums(nums)
    op_syms = symbol_op.(nums)
    #vec_of_pairs = Vector{Pair{Symbol, Union{Int,UnitRange{Int}}}}()
    vec_of_pairs = Vector{Pair{Symbol, UnitRange{Int}}}()
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

"""
    axis_of_nums(nums)

Return an axis with keys corresponding to each symbolic array and values
to the indices in the given vector of nums.
"""
function axis_of_nums(nums)
    pos_nums = indices_of_nums(nums)
    pos_nums_scalar = [first(pr) => length(last(pr)) == 1 ? first(last(pr)) : last(pr)
                       for pr in pos_nums]
    #CA._component_axis(first(axes(CA.ComponentArray(;pos_nums_scalar...))))
    first(getaxes(CA.ComponentArray(; pos_nums_scalar...)))
end
