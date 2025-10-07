import ComponentArrays as CA
#import ComponentArrays: _get_index_axis

function axis_length(ax::AbstractAxis) 
    lastindex(ax) - firstindex(ax) + 1
end
axis_length(::FlatAxis) = 0
axis_length(::UnitRange) = 0
axis_length(ax::ViewAxis) = axis_length(ax.ax)
axis_length(::Shaped1DAxis{(0,)}) = 0

# function _get_axis(x::AbstractArray) 
#     @info("Providing Parameters as Array was deprecated for performance?")
#     # depr?: need a full-fledged axis
#     Axis(Tuple(i for i in symbol_op.(x)))
# end
# function _get_axis(x::NTuple{N}) where N
#     Axis(Tuple(i for i in symbol_op.(x))::NTuple{N})
# end
function _get_axis(x::Tuple)
    Axis(Tuple(i for i in symbol_op.(x)))
end
_get_axis(x::ComponentVector) = first(getaxes(x))
_get_axis(x::AbstractAxis) = x
_get_axis(x::CA.CombinedAxis) = CA._component_axis(x)

# TODO move to ComponentArrays.jl
# type piracy: https://github.com/jonniedie/ComponentArrays.jl/issues/141
#@inline CA.getdata(x::ComponentVector) = getfield(x, :data)
# twutz 2311: remove CombinedAxis but only work with AbstractAxis (getaxes vs axes)
#attach_axis(x::AbstractVector, ax::CA.CombinedAxis) = ComponentArray(x, (CA._component_axis(ax),))
attach_axis(x::AbstractVector, ax::AbstractAxis) = ComponentArray(x, (ax,))
#attach_axis(x::ComponentVector, ax::AbstractAxis) = ComponentArray(getdata(x), (ax,))
function attach_axis(x::ComponentVector, ax::AbstractAxis)
    ComponentArray(getfield(x, :data), (ax,))
end
attach_x_axis(x::ComponentMatrix, ax::AbstractAxis) = ComponentArray(x, (ax, FlatAxis()))

# function _get_index_axis(cv::ComponentVector, ax::AbstractAxis)
#     first(getaxes(cv)) == ax && return(cv) # no need to reassamble
#     # extract subvectors and reassamble
#     keys_ax = keys(ax)
#     # k = keys_ax[1]
#     tmp = map(keys_ax) do k
#         cvs = cv[KeepIndex(k)]
#         axs = ax[k].ax
#         @show cvs, axs, length(cvs), CA.last_index(ax)
#         !(axs isa CA.NullAxis) && CA.last_index(ax) != length(cvs) && error(
#             "Expect axis extracting component `$k` to extract n=$(length(cvs)) elements, " * 
#             "but was $axs.")
#         cvk = _get_index_axis(cvs, axs)
#         ComponentVector(getdata(cvk), Viewax)
#     end
#     vcat(tmp...)
# end
# _get_index_axis(x, ax::CA.NullorFlatAxis) = x
# # in order to extract entire component, do not specify subaxes
# # e.g. (a=1) to match entire (a=(a1=1, a2=2))
# _get_index_axis(cv::ComponentVector, ax::CA.NullorFlatAxis) = cv # else method ambiguous

# function _set_index_axis!(cv::ComponentVector, s::ComponentVector)
#     ax = first(getaxes(s))
#     if first(getaxes(cv)) == ax 
#         cv .= s
#         return(cv)        
#     end
#     # setindex in each direct component
#     keys_ax = keys(ax)
#     for k in keys_ax
#         cvs = getproperty(cv, k)
#         ss = getproperty(s,k)
#         #@show cvs, ss
#         if cvs isa ComponentVector
#             _set_index_axis!(cvs, ss)
#         else
#             setproperty!(cv, k, ss)
#         end
#     end
#     cv
# end
# # function _set_index_axis!(cv::ComponentVector, sany)
# #     cv .= sany
# # end
# i_test = () -> begin
#     cvs = pset.is_state[i] ? u0l[k] : (pset.is_state[i] ? pl[k] : missing) 

#     cv0 = ComponentVector(a=(a1=1,a2=2,a3=3),b=20)
#     cv = copy(cv0)
#     cvs = ComponentVector(a=(a1=1,a3=3))
#     s = cv0[cvs]
#     s.a.a3 = 6
#     _set_index_axis!(cv, s)

#     cvs = ComponentVector(a=1)
#     s = cv0[cvs]
#     s.a = [10,20,30]
#     _set_index_axis!(cv, s)    
# end

# """
# assembles a new ComponentVector with some components replaced by corresponding
# components of s.
# The Base array used is Vector which is converted to the grand*parent type of cv,
# where `type_of(parent(x)) = type_of(x)`. Only base-types that support conversion
# from Vector without additional information are supported, e.g. StaticArray, but
# not AxisArray or NamedArray.
# Use `_set_index_axis!(cv, s)` to mutate an existing vector and keep its type
# """
# function _update_cv(cv::ComponentVector{T,A}, s::ComponentVector) where {T,A}
#     ax = first(getaxes(cv))
#     ax == first(getaxes(s)) && return s
#     keys_ax = keys(ax)
#     g = (ks ∈ keys_ax for ks in keys(first(getaxes(s)))) 
#     if !all(g) 
#         ks_missing = keys(first(getaxes(s)))[.! g]
#         error(
#             "updating keys $(ks_missing) not found in updated component vector with keys $ax")
#     end
#     tmp = map(keys_ax) do k
#         cvs = getproperty(cv, k)
#         !haskey(s, k) && return(cvs)
#         ss = getproperty(s,k)
#         #@show cvs, ss
#         if cvs isa ComponentVector
#             (ss isa ComponentVector) && return(_update_cv(cvs, ss))
#             (length(ss) != length(cvs)) && error(
#                 "Expected updating argument `s` to match ComponentVector $k=$cvs of " * 
#                 "length($(length(cvs))), but was $ss.")
#             attach_axis(ss, first(getaxes(cvs)))
#         else
#             (length(ss) != length(cvs)) && error(
#                 "Expected updating key $k of length($(length(cvs))), with argument $ss.")
#             ss
#         end
#     end
#     res = ComponentVector(NamedTuple{keys_ax}(tmp))
#     # make sure to use the same base-type instead of vector
#     # direct parent maybe SubArray, which one cannot convert to, need to call parent(parent(...))
#     AT = typeof(fixpoint(parent,cv;fmap=typeof))
#     #AT = typeof(parent(cv))
#     attach_axis(convert(AT,getdata(res))::AT, first(getaxes(cv)))
# end

# function _ax_string_prefixed(any; prefix);  ""; end
# _ax_string_prefixed(sym::Symbol; prefix) = [string(sym)]
# function _ax_string_prefixed(v::AbstractArray; prefix); "[" .* string.(eachindex(v)) .* "]" end
# function _ax_string_prefixed(ax::CA.ShapedAxis; prefix="₊") 
#     #@show "ShapedAxis", ax
#     if length(CA.indexmap(ax)) != 0
#         error("Unanticipated case: ax")
#     else
#         "[" .* _enum_dim(size(ax)) .* "]"
#     end
# end
# function _ax_string_prefixed(ax::AbstractAxis; prefix="₊") 
#     ax_inner = Axis(ax)
#     if ax_inner == ax
#         _ax_string_prefixed(CA.indexmap(ax); prefix)
#     else
#         _ax_string_prefixed(ax_inner; prefix)
#     end
# end
# # use getaxes rather than axes to extract AbstractAxis instead of CombinedAxis
# # function _ax_string_prefixed(ax::CA.CombinedAxis; prefix="₊") 
# #     #@show "combinedAxis"
# #     _ax_string_prefixed(CA._component_axis(ax); prefix)
# # end
# #_ax_symbol(t::Tuple{Vararg{<:AbstractAxis}}) = mapreduce(vcat, _ax_symbol, t)
# function _ax_string_prefixed(t::NamedTuple; prefix="₊") 
#     #@show "NamedTuple", t, prefix
#     tmp = [prefix .* string(k) .* _ax_string_prefixed(v; prefix) for (k,v) in pairs(t)]
#     reduce(vcat, tmp; init=[])
# end 
# #_extend = (s1,s2) -> isnothing(s2) ? s1 : s1 * "₊" * s2
# #_ax_symbols(ax) = Tuple(Symbol.(getindex.(_ax_string_prefixed(ax), Ref(2:end))))

# "Generate all cartesian indexing combinations"
# function _enum_dim(s,d) 
#     d == 1 && return(string.(1:s[1]))
#     reduce(vcat, (_enum_dim(s,d-1) .* "," .* string(i) for i in 1:s[d]))
# end
# function _enum_dim(s) 
#     _enum_dim(s, length(s))
# end


"""
    _ax_symbols_tuple(ax::AbstractAxis; prefix = "₊")

Expand Vector axes such as `ViewAxis(1:3, Axis(a = 1:3,))` 
to `Symbol.(("a[1]","a[2]","a[3]"))`
"""
function _ax_symbols_tuple(ax::AbstractAxis; prefix = "₊")
    (labels_noprefix(ax) .|> (x -> replace(x, "." => prefix)) .|> Symbol) |> Tuple
end
function _ax_symbols_tuple(ax::UnitRange) # representing a 0.length FlatAxis
    ()
end
function _ax_symbols_vector(ax::AbstractAxis; prefix = "₊")
    (labels_noprefix(ax) .|> (x -> replace(x, "." => prefix)) .|> Symbol)::Vector{Symbol}
end

# CA.labels(ax) would be type piracy II 
labels_noprefix(ax::AbstractAxis) = map(x -> x[2:end], _labels(ax))

function _labels(x::AbstractAxis, nview::Int = 0)
    #@info "Absract:$(typeof(x))"
    vcat((".$(key)" .* _labels(x[key]) for key in keys(x))...)
end
function _labels(x::NamedTuple, nview::Int = 0)
    length(x) == 0 ? [""] : vcat((".$(key)" .* _labels(x[key]) for key in keys(x))...)
end

function _labels(x::AbstractArray, nview::Int = 0)
    vcat(("[" * join(i.I, ",") * "]" for i in CartesianIndices(x))...)
end
function _labels(x, nview::Int = 0)
    # @info "_labels(x)"
    # @show x, typeof(x)
    ""
end
_labels(x::CA.NullAxis, nview::Int = 0) = ""

function _labels(x::CA.PartitionedAxis{PartSz, IdxMap}, nview::Int) where {PartSz, IdxMap}
    la = _labels(IdxMap)
    ncomp = Int(nview / PartSz)
    #@show nview, PartSz, ncomp
    v = vcat(("[$i]" for i in 1:ncomp)...)
    vcat((vi .* a for (a, vi) in Iterators.product(la, v))...)
end
#function _labels(x::CA.ShapedAxis{Shape, IdxMap}, nview::Int) where {Shape, IdxMap}
function _labels(x::CA.ShapedAxis{Shape}, nview::Int) where {Shape}
    #la = _labels(IdxMap)
    # idx = CA.indexmap(x)
    # @show typeof(idx), idx
    #la = _labels(CA.indexmap(x))
    vcat(("[" * join(i.I, ",") * "]" for i in CartesianIndices(Shape))...)
    #@show la, v
    #vcat((vi .* a for (a, vi) in Iterators.product(la, v))...)
end

function _labels(x::CA.ComponentIndex{N, FlatAxis}, nview::Int = 0) where {N}
    vcat(("[$i]" for i in eachindex(x.idx))...)
end
function _labels(x::CA.ComponentIndex{N, <:AbstractAxis}, nview::Int = 0) where {N}
    _labels(x.ax, length(x.idx))
end



"""
returns a plain new SVector.
"""
function getindex_svector(cv::ComponentVector, k::Symbol)
    cv_k = @view(cv[KeepIndex(k)])
    SVector{axis_length(first(getaxes(cv_k)))}(getdata(cv_k))
    # only available with ComponentArrays 0.15, need to downgrade to 0.13.8 to use Turing
    #getdata(CA.static_getproperty(cv, Val(k)))
end

function _update_cv_top(cv::ComponentVector, s::ComponentVector)
    keyss = keys(s)
    mkeys = (!(k ∈ keys(cv)) for k in keyss)
    any(mkeys) && error("The following keys to update were not found in destination " *
          string(keyss[collect(mkeys)]))
    gen_is_updated = (k ∈ keys(s) for k in keys(cv))
    #is_updated = collect(gen_is_updated)
    is_updated = SVector{length(keys(_get_axis(cv)))}(gen_is_updated...)
    _update_cv_top(cv, s, is_updated)
end

"""
    _update_cv_top(cv::ComponentVector{TD}, s::ComponentVector{TS}, is_updated)

Return a new ComponentVector of eltype `promote_type(TD, TS)` with those components 
in `cv` at position, `i`, for which `is_key_updated[i]` is true, are replaced by the 
corresponding name of source s. 
"""
function _update_cv_top(cv::ComponentVector{TD, TAD},
        s::ComponentVector{TS},
        is_updated::AbstractVector{Bool}) where {TD, TAD, TS}
    # s has not entries, return a copy
    axis_length(_get_axis(s)) == 0 && copy(cv)
    T_EL = promote_type(TD, TS)
    #(i,k) = first(enumerate(keys(cv)))
    #(i,k) = last(enumerate(keys(cv)))
    ftmp = (i, k) -> begin
        if is_updated[i]
            # extracting the underlying array does not gain performance but 
            # makes problems in vcat? 
            #val = @view s[k]
            #val_s = @view s[KeepIndex(k)]
            val_s = getindex_svector(s, k)
        else
            #val_cv = @view cv[KeepIndex(k)]
            val_cv = getindex_svector(cv, k)
        end
    end
    g = (ftmp(i, k) for (i, k) in enumerate(keys(cv)))
    #data = vcat(g...) # does not work with julia 1.6 - freezes
    data = reduce(vcat, g, init = StaticArrays.SA[])
    T_C = TAD <: StaticArray ?
          similar_type(TAD, T_EL) :
          typeof(similar(getdata(cv), T_EL))
    data_conv = convert(T_C, data)::T_C
    attach_axis(data_conv, _get_axis(cv))
end

# # TODO: wait for ComponentArrays implement length(Axis)
# """
#     subaxis(ax, sym::Symbol)
#     subaxis(ax, syms)

# Construct new reindexed axis for a top-level subcomponent of an axis.
# """
# function subaxis(ax::AbstractAxis, syms); _subaxis(ax, syms); end,
# function subaxis(ax::AbstractAxis, sym::Symbol); _subaxis(ax, (sym,)); end,
# function subaxis(cv::ComponentVector, syms); subaxis(first(getaxes(cv)), syms); end

# function _subaxis(ax::Axis,syms)
#     is_missing = map(s -> !(s ∈ keys(ax)), syms)
#     any(is_missing) && error(
#         "Expected subcomponents to be among keys(ax)=$(keys(ax)). Failed for " * 
#         "$([s for (m,s) in zip(is_missing, syms) if m])")
#     length_axs = NamedTuple{keys(CA.indexmap(ax))}(map(length, CA.indexmap(ax)))
#     # start positions of subaxes in original and in target axis
#     start_axs = NamedTuple{keys(length_axs)}(
#         cumsum(vcat(1,collect(length_axs)[1:end-1])))
#     start_axt = NamedTuple{syms}(cumsum(
#         vcat(1,[p.second for p in pairs(length_axs) if p.first ∈ syms])[1:end-1]))
#     nts = map(syms) do sym
#         ax_sym = indexmap(ax)[sym]
#         reindex(ax_sym, start_axt[sym]-start_axs[sym])
#     end"Y[1]"
#     Axis(NamedTuple{syms}(nts))
# end

"""
Removes the highest level of keys.
Keeps the reference to the underlying data, but changes the axis.
If first-level vector has no sub-names, an error (Aguement Error tuple must be non-empty)
is thrown.
"""
function flatten1(cv::ComponentVector)
    # return a tuple of (key, value) as zip(keys, values) would do
    # gen = (((ks, cv[k][ks]) for ks in keys(cv[k])) for k in keys(cv) if !isempty(cv[k]))
    # cv_new = CA.ComponentVector(; Iterators.Flatten(gen)...)
    # benchmarks show that the vcat variant is more efficient
    if length(cv)==0 
        return attach_axis(cv, FlatAxis())
    else
        gen_cvs = (cv[k] for k in keys(cv) if !isempty(cv[k]))
        cv_new = reduce(vcat, gen_cvs)
        attach_axis(cv, first(getaxes(cv_new)))
    end
end

"""
Takes several ComponentVectors of subvectors and merges each subvector
MAYBE
"""
function merge_subvectors(cvs...;mkeys)
    # k = :state # k = :par
    # allkeys = union(collect.(keys.(cvs))...) # unclear order
    tup = map(mkeys) do k
        tmp = (x[k] for x in cvs if k ∈ keys(x))
        tmp2 = reduce(vcat, tmp)
        CA.ComponentVector(;zip((k,), (tmp2,))...)
    end
    tup_nonempty = (t for t in tup if !isempty(t))
    reduce(vcat,tup_nonempty)
end
#popt = merge_subvectors(fixed, random, indiv; mkeys=(:state, :par))

"""
    map_keys(FUN, cv::ComponentVector; rewrap::Val{is_rewrap}=Val(true)) 

Apply a function to all sub-vectors and create a new ComponentVector with
the same top-level keys.
"""
function map_keys(FUN, cv::ComponentVector; rewrap::Val{is_rewrap}=Val(true)) where is_rewrap
    tup = map(keys(cv)) do k
        ret = FUN(cv[k])
        eltype(ret) != Union{} || error("For mapping empty keys, provide a proper eltype."*
        "Did you accidentally write key=[] instead of e.g. Float64[] ?")
        ret
    end
    is_rewrap ? ComponentVector(;zip(keys(cv),tup)...) : tup
end
