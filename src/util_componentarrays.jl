import ComponentArrays as CA

function _get_index_axis(cv::ComponentVector, ax::AbstractAxis)
    first(getaxes(cv)) == ax && return(cv) # no need to reassamble
    # extract subvectors and reassamble
    keys_ax = keys(ax)
    # k = keys_ax[1]
    tmp = map(keys_ax) do k
        cvs = getproperty(cv, k)
        axs = ax[k].ax
        #@show cvs, axs, length(cvs), CA.last_index(ax)
        !(axs isa CA.NullAxis) && CA.last_index(ax) != length(cvs) && error(
            "Expect axis extracting component `$k` to extract n=$(length(cvs)) elements, " * 
            "but was $axs.")
        _get_index_axis(cvs, axs)
    end
    ComponentVector(NamedTuple{keys_ax}(tmp))
end
_get_index_axis(x, ax::CA.NullorFlatAxis) = x
# in order to extract entire component, do not specify subaxes
# e.g. (a=1) to match entire (a=(a1=1, a2=2))
_get_index_axis(cv::ComponentVector, ax::CA.NullorFlatAxis) = cv # else method ambiguous

function _set_index_axis!(cv::ComponentVector, s::ComponentVector)
    ax = first(getaxes(s))
    if first(getaxes(cv)) == ax 
        cv .= s
        return(cv)        
    end
    # setindex in each direct component
    keys_ax = keys(ax)
    for k in keys_ax
        cvs = getproperty(cv, k)
        ss = getproperty(s,k)
        #@show cvs, ss
        if cvs isa ComponentVector
            _set_index_axis!(cvs, ss)
        else
            setproperty!(cv, k, ss)
        end
    end
    cv
end
# function _set_index_axis!(cv::ComponentVector, sany)
#     cv .= sany
# end
i_test = () -> begin
    cvs = pset.is_state[i] ? u0l[k] : (pset.is_state[i] ? pl[k] : missing) 

    cv0 = ComponentVector(a=(a1=1,a2=2,a3=3),b=20)
    cv = copy(cv0)
    cvs = ComponentVector(a=(a1=1,a3=3))
    s = cv0[cvs]
    s.a.a3 = 6
    _set_index_axis!(cv, s)


    cvs = ComponentVector(a=1)
    s = cv0[cvs]
    s.a = [10,20,30]
    _set_index_axis!(cv, s)    
end


"""
assembles a new ComponentVector with some components replaced by corresponding
components of s.
The Base array used is Vector which is converted to the grand*parent type of cv,
where `type_of(parent(x)) = type_of(x)`. Only base-types that support conversion
from Vector without additional information are supported, e.g. StaticArray, but
not AxisArray or NamedArray.
Use `_set_index_axis!(cv, s)` to mutate an exsiting vector and keep its type
"""
function _update_cv(cv::ComponentVector{T,A}, s::ComponentVector) where {T,A}
    ax = first(getaxes(cv))
    ax == first(getaxes(s)) && return s
    keys_ax = keys(ax)
    g = (ks ∈ keys_ax for ks in keys(first(getaxes(s)))) 
    if !all(g) 
        ks_missing = keys(first(getaxes(s)))[.! g]
        error(
            "updating keys $(ks_missing) not found in updated compoenent vector with keys $ax")
    end
    tmp = map(keys_ax) do k
        cvs = getproperty(cv, k)
        !haskey(s, k) && return(cvs)
        ss = getproperty(s,k)
        #@show cvs, ss
        if cvs isa ComponentVector
            (ss isa ComponentVector) && return(_update_cv(cvs, ss))
            (length(ss) != length(cvs)) && error(
                "Expected updating argument `s` to match ComponentVector $k=$cvs of " * 
                "length($(length(cvs))), but was $ss.")
            attach_axis(ss, first(getaxes(cvs)))
        else
            (length(ss) != length(cvs)) && error(
                "Expected updating key $k of length($(length(cvs))), with argument $ss.")
            ss
        end
    end
    res = ComponentVector(NamedTuple{keys_ax}(tmp))
    # make sure to use the same base-type instead of vector
    # direct parent maybe SubArray, which one cannot convert to, need to call parent(parent(...))
    AT = typeof(fixpoint(parent,cv;fmap=typeof))
    #AT = typeof(parent(cv))
    attach_axis(convert(AT,getdata(res))::AT, first(getaxes(cv)))
end


function _ax_string_prefixed(any; prefix);  ""; end
_ax_string_prefixed(sym::Symbol; prefix) = [string(sym)]
function _ax_string_prefixed(v::AbstractArray; prefix); "[" .* string.(eachindex(v)) .* "]" end
function _ax_string_prefixed(ax::CA.ShapedAxis; prefix="₊") 
    #@show "ShapedAxis", ax
    if length(CA.indexmap(ax)) != 0
        error("Unanticipated case: ax")
    else
        "[" .* _enum_dim(size(ax)) .* "]"
    end
end
function _ax_string_prefixed(ax::AbstractAxis; prefix="₊") 
    ax_inner = Axis(ax)
    if ax_inner == ax
        _ax_string_prefixed(CA.indexmap(ax); prefix)
    else
        _ax_string_prefixed(ax_inner; prefix)
    end
end
# use getaxes rather than axes to extract AbstractAxis instead of CombinedAxis
# function _ax_string_prefixed(ax::CA.CombinedAxis; prefix="₊") 
#     #@show "combinedAxis"
#     _ax_string_prefixed(CA._component_axis(ax); prefix)
# end
#_ax_symbol(t::Tuple{Vararg{<:AbstractAxis}}) = mapreduce(vcat, _ax_symbol, t)
function _ax_string_prefixed(t::NamedTuple; prefix="₊") 
    #@show "NamedTuple", t, prefix
    tmp = [prefix .* string(k) .* _ax_string_prefixed(v; prefix) for (k,v) in pairs(t)]
    reduce(vcat, tmp; init=[])
end 
#_extend = (s1,s2) -> isnothing(s2) ? s1 : s1 * "₊" * s2
#_ax_symbols(ax) = Tuple(Symbol.(getindex.(_ax_string_prefixed(ax), Ref(2:end))))

"Generate all cartesian indexing combinations"
function _enum_dim(s,d) 
    d == 1 && return(string.(1:s[1]))
    reduce(vcat, (_enum_dim(s,d-1) .* "," .* string(i) for i in 1:s[d]))
end
function _enum_dim(s) 
    _enum_dim(s, length(s))
end



