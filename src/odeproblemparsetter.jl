"""
    ODEProblemParSetter(state_template,par_template,popt_template) 
    ODEProblemParSetter(sys::ODESystem, popt_template) 

Helps keeping track of a subset of initial states and parameters to be optimized.

# Arguments
- `state_template`: ComponentVector or Axis of all the initial states of the problem
- `par_template`: all the parameters of the problem
- `popt_template`: the parameter/initial states to be optimized. 
  If given as Tuple or AbstractVector of symbols, then a template ComponentVector
  is extracted from `state_template` and `par_template`.

If all of `state_template`, `par_template`, and `popt_template` are type-inferred Axes,
then also the constructed ODEProblemParSetter is type-inferred.

The states and parameters can be extracted from an `ModelingToolkit.ODESystem`.

Note the similar [`ODEProblemParSetterConcrete`](@ref) with template parameters, which 
supports type-stable calls (see [Concrete ProblemUpdater](@ref)).
"""
struct ODEProblemParSetter <: AbstractODEProblemParSetter
    ax_paropt::AbstractAxis
    ax_state::AbstractAxis
    ax_par::AbstractAxis
    is_updated_state_i::AbstractVector
    is_updated_par_i::AbstractVector
    ax_paropt_scalar::AbstractAxis
    ax_paropt_flat1::AbstractAxis
    function ODEProblemParSetter(ax_state::AbstractAxis,
            ax_par::AbstractAxis, ax_paropt::AbstractAxis, ax_paropt_scalar::AbstractAxis,
            ax_paropt_flat1::AbstractAxis, ::Val{isval}) where {isval}
        if isval
            is_valid, msg = validate_keys_state_par(ax_paropt_scalar, ax_state, ax_par)
            !is_valid && error(msg)
        end
        keys_paropt_state = keys(CA.indexmap(ax_paropt_scalar)[:state])
        keys_paropt_par = keys(CA.indexmap(ax_paropt)[:par])
        is_updated_state_i = isempty(keys_paropt_state) ?
                             SVector{0, Bool}() :
                             SVector((k ∈ keys_paropt_state for k in keys(ax_state))...)
        is_updated_par_i = isempty(keys_paropt_par) ?
                           SVector{0, Bool}() :
                           SVector((k ∈ keys_paropt_par for k in keys(ax_par))...)
        new(ax_paropt, ax_state, ax_par, is_updated_state_i, is_updated_par_i,
            ax_paropt_scalar, ax_paropt_flat1)
    end
end

function ODEProblemParSetter(state_template, par_template, popt_template,
        system::Union{AbstractODESystem, Nothing} = nothing;
        is_validating = Val{true}())
    ax_paropt = _get_axis(popt_template)
    ax_state = _get_axis(state_template)
    ax_par = _get_axis(par_template)
    ax_state_array = isnothing(system) ? ax_state : axis_of_nums(states(system))
    if !(:state ∈ keys(ax_paropt) || :par ∈ keys(ax_paropt))
        ax_paropt = assign_state_par(ax_state_array, ax_par, ax_paropt)
    end
    (ax_state_scalar, ax_paropt_scalar) = isnothing(system) ?
                                          (ax_state, ax_paropt) :
                                          scalarize_par_and_paroptstate(ax_state,
        ax_paropt, system)
    #Main.@infiltrate_main        
    ax_paropt_flat1 = try
        first(getaxes(flatten1(ComponentVector(1:axis_length(ax_paropt),ax_paropt))))
    catch e
        @warn("Could not produce flattened version of axis_paropt=$ax_paropt. " * 
        "Are there duplicate keys in state and parameters?")
        FlatAxis()
    end
    ODEProblemParSetter(ax_state_scalar, ax_par, ax_paropt, ax_paropt_scalar, ax_paropt_flat1, is_validating)
end

"scalarize vector-valued entries in state and paropt.state"
function scalarize_par_and_paroptstate(ax_state::AbstractAxis,
        ax_paropt::AbstractAxis,
        system::AbstractODESystem)
    cv = ComponentArray(1:axis_length(ax_state), ax_state)
    cvs = expand_base_num_axes(cv, system)
    ax_state_scalar = first(getaxes(cvs))
    cv = ComponentArray(1:axis_length(ax_paropt), ax_paropt)
    cvs = ComponentVector(state = expand_base_num_axes(cv.state, system), par = cv.par)
    ax_paropt_scalar = first(getaxes(cvs))
    (ax_state_scalar, ax_paropt_scalar)
end

function ODEProblemParSetter(state_template,
        par_template, popt_template::Union{NTuple{N, Symbol}, AbstractVector{Symbol}},
        system::Union{AbstractODESystem, Nothing} = nothing;
        is_validating = Val{true}()) where {N}
    ax_par = _get_axis(par_template)
    ax_state_array = isnothing(system) ?
                     _get_axis(state_template) :
                     axis_of_nums(states(system))
    popt_template_new = length(popt_template) == 0 ? Axis() : begin
        # construct a template by extracting the components of u0 and p
        u0 = attach_axis(1:axis_length(ax_state_array), ax_state_array)
        p = attach_axis(1:axis_length(ax_par), ax_par)
        u0p = vcat(u0, p)
        u0p isa ComponentArray || error("Could not concatenate u0=$u0 and p=$p.")
        popt_template_new = u0p[popt_template]
    end
    ODEProblemParSetter(state_template, ax_par, popt_template_new, system; is_validating)
end

function assign_state_par(ax_state, ax_par, ax_paropt)
    state_keys = Vector{Symbol}()
    par_keys = Vector{Symbol}()
    for key in keys(ax_paropt)
        key ∈ keys(ax_state) && push!(state_keys, key)
        key ∈ keys(ax_par) && push!(par_keys, key)
    end
    missing_keys = setdiff(keys(ax_paropt), vcat(state_keys, par_keys))
    length(missing_keys) != 0 && @warn("Expected optimization parameters to be part of "*
          "state or parameters, but did not found parameters "*string(missing_keys)*".")
    duplicate_keys = intersect(state_keys, par_keys)
    length(duplicate_keys) != 0 && @warn("Expected optimization parameters to be either "*
          " part of state or parameters, but following occur in both "*
          string(duplicate_keys)*". Will update those only in state.")
    # assume to refer to state only
    par_keys = setdiff(par_keys, duplicate_keys)
    tmp = attach_axis((1:axis_length(ax_paropt)), ax_paropt)
    # empty ComponentVector does not translate to tmp2
    # tmp_state = isempty(state_keys) ? ComponentVector() : @view tmp[state_keys]
    # tmp_par = isempty(par_keys) ? ComponentVector() : @view tmp[par_keys]
    tmp_state = @view tmp[state_keys]
    tmp_par = @view tmp[par_keys]
    tmp2 = CA.ComponentVector(state = tmp_state, par = tmp_par)
    (state_keys..., par_keys...) != keys(ax_paropt) &&
        @warn("expected ax_paropt to contain state keys first: "*
              "$((state_keys..., par_keys...)), "*"but was $(keys(ax_paropt)). " *
              "label_paropt(), get_paropt(), etc. assume the first order.")
    return _get_axis(tmp2)
end

function ODEProblemParSetter(sys::ODESystem, paropt)
    # ODEProblemParSetter(axis_of_nums(states(sys)), axis_of_nums(parameters(sys)), paropt, sys)
    #ODEProblemParSetter(Axis(Symbol.(states(sys))), axis_of_nums(parameters(sys)), paropt, sys)
    # simplify X(t) to X but keep (Y(t))[i] intact
    ODEProblemParSetter(Axis(symbol_op_scalar.(states(sys))),
        axis_of_nums(parameters(sys)),
        paropt,
        sys)
end

function get_concrete(pset::ODEProblemParSetter)
    ODEProblemParSetterConcrete(pset.ax_state,
        pset.ax_par,
        pset.ax_paropt,
        pset.ax_paropt_scalar,
        pset.ax_paropt_flat1,
        Val{false}())
end

ODEProblemParSetterU = Union{ODEProblemParSetter, ODEProblemParSetterConcrete}

axis_state(ps::ODEProblemParSetterU) = ps.ax_state
axis_par(ps::ODEProblemParSetterU) = ps.ax_par
axis_paropt(ps::ODEProblemParSetterU) = ps.ax_paropt
axis_paropt_scalar(ps::ODEProblemParSetterU) = ps.ax_paropt_scalar
axis_paropt_flat1(ps::ODEProblemParSetterU) = ps.ax_paropt_flat1

classes_paropt(pset::ODEProblemParSetterU) = (:state, :par)

# # Using unexported interface of ComponentArrays.axis, one place to change
# "Accessor function for index from ComponentIndex"
# idx(ci::CA.ComponentIndex) = ci.idx

"""
    u0new, pnew = update_statepar(pset::ODEProblemParSetter, popt, u0, p) 

Return an updated problem or updates states and parameters where
values corresponding to positions in `popt` are set.
"""
function update_statepar(pset::ODEProblemParSetterU, popt, u0, p)
    poptc = attach_axis(popt, axis_paropt_scalar(pset))
    u0c = attach_axis(u0, axis_state(pset)) # no need to copy 
    pc = attach_axis(p, axis_par(pset))
    u0_new = isempty(poptc.state) ?
             copy(u0c) :
             _update_cv_top(u0c, poptc.state, pset.is_updated_state_i)
    p_new = isempty(poptc.par) ?
            copy(pc) :
            _update_cv_top(pc, poptc.par, pset.is_updated_par_i)
    # if u0 was not a ComponentVector, return then data inside
    return (u0 isa ComponentVector) ? u0_new : getfield(u0_new, :data),
    (p isa ComponentVector) ? p_new : getfield(p_new, :data)
end

# mutating version does not work with derivatives
# function update_statepar(pset::ODEProblemParSetter, popt, u0, p)
#     poptc = attach_axis(popt, axis_paropt(pset))
#     u0c = attach_axis(copy(u0), axis_state(pset))
#     pc = attach_axis(copy(p), axis_par(pset))
#     for k in keys(poptc.state)
#         u0c[k] = poptc.state[k]
#     end
#     for k in keys(poptc.par)
#         pc[k] = poptc.par[k]
#     end
#     # if u0 was not a ComponentVector, return then data inside
#     return (u0 isa ComponentVector) ? u0c : getfield(u0c, :data), 
#         (p isa ComponentVector) ? pc : getfield(pc, :data)
# end

# struct ODEVectorCreator <: AbstractVectorCreator; end
# function (vc::ODEVectorCreator)(pset, u0, p)
#     ET = promote_type(eltype(u0), eltype(p))
#     Vector{ET}(undef, count_paropt(pset))
# end

# struct ODEMVectorCreator <: AbstractVectorCreator; end
# function (vc::ODEMVectorCreator)(pset, u0, p)
#     ET = promote_type(eltype(u0), eltype(p))
#     N = count_paropt(pset)
#     MVector{N,ET}(undef)
# end

"""
    get_paropt_labeled(pset::ODEProblemParSetter, u0, p)

Returns a ComponentVector filled with corresponding entries from u0 and p.

Underlying type corresponds defaults to Vector. 
It cannot be fully inferred from u0 and p, because their type may hold
additional structure, such as length or names.
For obtaining a StaticVector instead, pass MTKHelpers.ODEMVectorCreator()
as the fourth argument. This method should work with any AbstractVectorCreator
that returns a mutable AbstractVector.
"""
function get_paropt_labeled(pset::ODEProblemParSetterU, u0, p;
        flat1::Val{is_flatten1} = Val(false)
        #vec_creator::AbstractVectorCreator=ODEVectorCreator()
) where is_flatten1
    u0c = attach_axis(u0, axis_state(pset))
    pc = attach_axis(p, axis_par(pset))
    ax = axis_paropt_scalar(pset)
    k_state = keys(CA.indexmap(ax).state)
    k_par = keys(CA.indexmap(ax).par)
    #
    # axis is attached later - use SVector instead of view into CA to avoid recompiling vcat
    gen_state = (getindex_svector(u0c, k) for k in k_state)
    gen_par = (getindex_svector(pc, k) for k in k_par)
    # gen_state = (@view(u0c[KeepIndex(k)]) for k in k_state)
    # gen_par = (@view(pc[KeepIndex(k)]) for k in k_par)
    #
    # gen_state = (u0c[k] for k in k_state) # more allocations
    # gen_par = (pc[k] for k in k_par)
    # gen_state = (@view(u0c[k]) for k in k_state) # takes long than with KeepIndex
    # gen_par = (@view(pc[k]) for k in k_par)
    # Main.@infiltrate_main
    # tmp = collect(gen_state)
    # tmp = collect(gen_par)
    #_data = vcat(gen_state..., gen_par...) # does not work with Julia 1.6
    _data = vcat(reduce(vcat, gen_state, init = StaticArrays.SA[]),
        reduce(vcat, gen_par, init = StaticArrays.SA[]))
    T = promote_type(eltype(u0), eltype(p))
    paropt = (
        is_flatten1 ? label_paropt_flat1(pset, _data) : label_paropt(pset, _data)
        )::ComponentVector{T, AV} where {AV <: AbstractVector{T}}
    # cv_state = ComponentVector((;zip(k_state,(u0c[k] for k in k_state))...))
    # cv_par = ComponentVector((;zip(k_par,(pc[k] for k in k_par))...))
    # paropt = ComponentVector(state = cv_state, par = cv_par)
    return paropt
end

# mutating version does not work with gradient
# function get_paropt_labeled(pset::ODEProblemParSetter, u0, p, 
#     vec_creator::AbstractVectorCreator=ODEVectorCreator())
#     u0c = attach_axis(u0, axis_state(pset))
#     pc = attach_axis(p, axis_par(pset))
#     #Main.@infiltrate_main
#     data = vec_creator(pset, u0, p)
#     paropt = attach_axis(data, axis_paropt(pset))
#     for k in keys(paropt.state)
#         paropt.state[k] = u0c[k]
#     end
#     for k in keys(paropt.par)
#         paropt.par[k] = pc[k]
#     end
#     return paropt
# end

"""
    validate:keys(pset)

Checks whether all components of paropt-Axis are occurring
in corresponding axes.     
"""
function validate_keys(pset::ODEProblemParSetterU)
    validate_keys_state_par(axis_paropt_scalar(pset), axis_state(pset), axis_par(pset))
end

function validate_keys_state_par(ax_paropt::AbstractAxis,
        ax_state::AbstractAxis,
        ax_par::AbstractAxis)
    paropt = attach_axis((1:axis_length(ax_paropt)), ax_paropt)
    if keys(ax_paropt) != (:state, :par)
        return (; isvalid = false,
            msg = String127("Expected paropt to have classification keys (:state,:par) " *
                            "but was " * string(keys(paropt))))
    end
    paropt.state isa CA.ComponentVector ||
        length(paropt.state) == 0 ||  # special case of empty ComponentVector, e.g. 2:1
        return return (; isvalid = false,
            msg = String127("Expected paropt.state <: ComponentVector, but was not."))
    paropt.par isa CA.ComponentVector ||
        length(paropt.par) == 0 ||  # special case of empty ComponentVector
        return return (; isvalid = false,
            msg = String127("Expected paropt.par <: ComponentVector, but was not."))
    u0c = attach_axis((1:axis_length(ax_state)), ax_state)
    for k in keys(paropt.state)
        k ∉ keys(u0c) && return (; isvalid = false,
            msg = String127("Expected optimined paropt.state.$k to be part of state, " *
                            "but was not."))
        length(paropt.state[k]) != length(u0c[k]) && return (; isvalid = false,
            msg = String127("Expected optimized paropt.state.$k to be of length " *
                            "$(length(u0c[k])) but had length $(length(paropt.state[k]))"))
    end
    pc = attach_axis((1:axis_length(ax_par)), ax_par)
    for k in keys(paropt.par)
        length_par_k = length(paropt.par[k])::Int
        k ∉ keys(pc) && return (; isvalid = false,
            msg = String127("Expected optimined paropt.par.$k to be part of parameters, " *
                            "but was not."))
        length_par_k != length(pc[k]) && return (; isvalid = false,
            msg = String127("Expected optimized paropt.par.$k to be of length " *
                            "$(length(pc[k])) but had length $length_par_k"))
    end
    return (; isvalid = true, msg = String127(""))
end

"""
    vcat_statesfirst(cvs; system)

Concatenate ComponentVectors, cvs, but move state entries before parameter entries.   
"""
function vcat_statesfirst(cvs...; system::AbstractSystem)
    popt0 = reduce(vcat,cvs)
    pset = LoggingExtras.withlevel(Logging.Error) do
        pset = ODEProblemParSetter(system, popt0)
    end
    popt_template = label_paropt_flat1(pset, 1:count_paropt(pset)) 
    popt0[keys(popt_template)]
end

