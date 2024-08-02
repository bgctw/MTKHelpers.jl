# type alias to save typing
# already defined in odeproblemparsetterconcrete.jl
#const VN = AbstractVector{<:SymbolicUtils.BasicSymbolic}

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
    ax_state_scalar::AbstractAxis
    ax_par::AbstractAxis
    ax_paropt_scalar::AbstractAxis
    ax_paropt_flat1::AbstractAxis
    opt_state_nums::VN
    opt_par_nums::VN
    par_ind::AbstractVector # propb.ps -> vector
    stateopt_ind::AbstractVector # to directly index into prob.u0
    popt_ind::AbstractVector # to directly index into prop.ps
    #state_scalars::AbstractDict{Symbol, AbstractVector}
    function ODEProblemParSetter(ax_state::AbstractAxis, ax_state_scalar::AbstractAxis,
        ax_par::AbstractAxis, ax_paropt::AbstractAxis, ax_paropt_scalar::AbstractAxis,
        ax_paropt_flat1::AbstractAxis,
        opt_state_nums::VN, opt_par_nums::VN,
        par_ind, stateopt_ind::AbstractVector, popt_ind::AbstractVector,
    ) 
        # VN <: AbstractVector{<:SymbolicUtils.BasicSymbolic} || error("expected VN <: AbstractVector{<:SymbolicUtils.BasicSymbolic}, but was $VN")
        keys_paropt_state = keys(CA.indexmap(ax_paropt_scalar)[:state])
        keys_paropt_par = keys(CA.indexmap(ax_paropt)[:par])
        new(ax_paropt, ax_state, ax_state_scalar, ax_par,
            ax_paropt_scalar, ax_paropt_flat1, 
            opt_state_nums, opt_par_nums,
            par_ind, stateopt_ind, popt_ind, 
            )
    end
end

function ODEProblemParSetter(state_template, par_template, popt_template,
    system::AbstractODESystem;
    is_validating::Val{isval}=Val{true}()) where {isval}
    ax_paropt = _get_axis(popt_template)
    ax_state = _get_axis(state_template)
    ax_par = _get_axis(par_template)
    ax_state_array = axis_of_nums(unknowns(system))
    # for each unique state or parameter (maybe vector), 
    if !(:state ∈ keys(ax_paropt) || :par ∈ keys(ax_paropt))
        ax_paropt = assign_state_par(ax_state_array, ax_par, ax_paropt)
    end
    scalar_num_map = get_scalar_num_map(system)
    scalar_num_map_sym = Dict(symbol_op(k) => v for (k,v) in scalar_num_map)
    scalar_num_map_symsym = Dict(k => symbol_num_getindex.(v) for (k,v) in scalar_num_map_sym)
    (ax_state_scalar, ax_paropt_scalar) = scalarize_par_and_paroptstate(ax_state,
        ax_paropt, scalar_num_map)
    ax_paropt_flat1 = try
        first(getaxes(flatten1(ComponentVector(1:axis_length(ax_paropt), ax_paropt))))
    catch e
        @warn("Could not produce flattened version of axis_paropt=$ax_paropt. " *
              "Are there duplicate keys in state and parameters?")
        FlatAxis()
    end
    if isval
        is_valid, msg = validate_keys_state_par(ax_paropt_scalar, ax_state_scalar, ax_par)
        !is_valid && error(msg)
    end
    _dict_nums = get_system_symbol_dict(system)
    par_nums = axis_length(ax_par) == 0 ? eltype(values(_dict_nums))[] : vcat(
        Symbolics.scalarize.(
            getindex.(Ref(_dict_nums), keys(ax_par)))...)
    T = eltype(par_nums)
    sym_paropt_par = keys(ax_paropt[:par].ax)
    opt_par_nums = (length(sym_paropt_par) == 0 ?
                    T[] :
                    vcat(
        getindex.(Ref(scalar_num_map_sym), sym_paropt_par)...))::Vector{T}
    sym_paropt_state = keys(ax_paropt[:state].ax)
    opt_state_nums = length(sym_paropt_state) == 0 ?
                     T[] :
                     vcat(
        getindex.(Ref(scalar_num_map_sym), sym_paropt_state)...)::Vector{T}
    par_ind = SII.parameter_index.(Ref(system), par_nums)
    stateopt_ind = SII.variable_index.(Ref(system), opt_state_nums)
    popt_ind = SII.parameter_index.(Ref(system), opt_par_nums)
    # provide ax_state_scalar instead of ax_state to match it to ax_paropt.state
    #   where only a subset of the vector indices can be referenced
    ODEProblemParSetter(
        ax_state, ax_state_scalar,
        ax_par, ax_paropt, ax_paropt_scalar, ax_paropt_flat1,
        opt_state_nums, opt_par_nums,
        par_ind, stateopt_ind, popt_ind, 
        )
end

"scalarize vector-valued entries in state and paropt.state"
function scalarize_par_and_paroptstate(ax_state::AbstractAxis,
    ax_paropt::AbstractAxis,
    scalar_num_map_sym::Dict{<:SymbolicUtils.BasicSymbolic,<:AbstractVector})
    scalar_num_map = Dict(symbol_op(k) => v for (k,v) in scalar_num_map_sym)
    cv = ComponentArray(1:axis_length(ax_state), ax_state)
    cvs = expand_base_num_axes(cv, scalar_num_map)
    ax_state_scalar = first(getaxes(cvs))
    cv = ComponentArray(1:axis_length(ax_paropt), ax_paropt)
    cvs = ComponentVector(state=expand_base_num_axes(cv.state, scalar_num_map), par=cv.par)
    ax_paropt_scalar = first(getaxes(cvs))
    (ax_state_scalar, ax_paropt_scalar)
end


function ODEProblemParSetter(state_template,
    par_template, popt_template::Union{NTuple{N,Symbol},AbstractVector{Symbol}},
    system::AbstractODESystem;
    is_validating=Val{true}()) where {N}
    ax_par = _get_axis(par_template)
    ax_state_array = isnothing(system) ?
                     _get_axis(state_template) :
                     axis_of_nums(unknowns(system))
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
    length(missing_keys) != 0 && @warn("Expected optimization parameters to be part of " *
                                       "state or parameters, but did not found parameters " * string(missing_keys) * ".")
    duplicate_keys = intersect(state_keys, par_keys)
    length(duplicate_keys) != 0 && @warn("Expected optimization parameters to be either " *
                                         " part of state or parameters, but following occur in both " *
                                         string(duplicate_keys) * ". Will update those only in state.")
    # assume to refer to state only
    par_keys = setdiff(par_keys, duplicate_keys)
    tmp = attach_axis((1:axis_length(ax_paropt)), ax_paropt)
    # empty ComponentVector does not translate to tmp2
    # tmp_state = isempty(state_keys) ? ComponentVector() : @view tmp[state_keys]
    # tmp_par = isempty(par_keys) ? ComponentVector() : @view tmp[par_keys]
    tmp_state = @view tmp[state_keys]
    tmp_par = @view tmp[par_keys]
    tmp2 = CA.ComponentVector(state=tmp_state, par=tmp_par)
    (state_keys..., par_keys...) != keys(ax_paropt) &&
        @warn("expected ax_paropt to contain state keys first: " *
              "$((state_keys..., par_keys...)), " * "but was $(keys(ax_paropt)). " *
              "label_paropt(), get_paropt(), etc. assume the first order.")
    return _get_axis(tmp2)
end

function ODEProblemParSetter(sys::ODESystem, paropt; is_validating=Val{true}())
    # ODEProblemParSetter(axis_of_nums(unknowns(sys)), axis_of_nums(parameters(sys)), paropt, sys)
    #ODEProblemParSetter(Axis(Symbol.(unknowns(sys))), axis_of_nums(parameters(sys)), paropt, sys)
    # simplify X(t) to X but keep (Y(t))[i] intact
    scalar_num_map = get_scalar_num_map(sys)
    ODEProblemParSetter(
        #Axis(symbol_op_scalar.(unknowns(sys))),
        axis_of_nums(unknowns(sys)),
        axis_of_nums(parameters(sys)),
        paropt,
        sys; is_validating)
end

function get_concrete(pset::ODEProblemParSetter)
    ODEProblemParSetterConcrete(pset.ax_state,
        pset.ax_state_scalar,
        pset.ax_par,
        pset.ax_paropt,
        pset.ax_paropt_scalar,
        pset.ax_paropt_flat1,
        pset.opt_state_nums, pset.opt_par_nums,
        pset.par_ind, pset.stateopt_ind, pset.popt_ind,   
    )
end

ODEProblemParSetterU = Union{ODEProblemParSetter,ODEProblemParSetterConcrete}

axis_state(ps::ODEProblemParSetterU) = ps.ax_state
axis_state_scalar(ps::ODEProblemParSetterU) = ps.ax_state_scalar
axis_par(ps::ODEProblemParSetterU) = ps.ax_par
axis_paropt(ps::ODEProblemParSetterU) = ps.ax_paropt
axis_paropt_scalar(ps::ODEProblemParSetterU) = ps.ax_paropt_scalar
axis_paropt_flat1(ps::ODEProblemParSetterU) = ps.ax_paropt_flat1

classes_paropt(::ODEProblemParSetterU) = (:state, :par)

function get_paropt(pset::ODEProblemParSetterU, prob::AbstractODEProblem; kwargs...) 
    u0_opt = prob.u0[pset.stateopt_ind]
    p_opt = prob.ps[pset.popt_ind]
    vcat(u0_opt, p_opt)
end

function get_par(pset::ODEProblemParSetterU, prob::SciMLBase.AbstractSciMLProblem)
    prob.ps[pset.par_ind]
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
    #validate_keys_state_par(axis_paropt(pset), axis_state(pset), axis_par(pset))
    validate_keys_state_par(axis_paropt_scalar(pset), axis_state_scalar(pset), axis_par(pset))
end

function validate_keys_state_par(ax_paropt::AbstractAxis,
    ax_state_scalar::AbstractAxis,
    ax_par::AbstractAxis)
    paropt = attach_axis((1:axis_length(ax_paropt)), ax_paropt)
    if keys(ax_paropt) != (:state, :par)
        return (; isvalid=false,
            msg=String127("Expected paropt to have classification keys (:state,:par) " *
                          "but was " * string(keys(paropt))))
    end
    paropt.state isa CA.ComponentVector ||
        length(paropt.state) == 0 ||  # special case of empty ComponentVector, e.g. 2:1
        return return (; isvalid=false,
            msg=String127("Expected paropt.state <: ComponentVector, but was not."))
    paropt.par isa CA.ComponentVector ||
        length(paropt.par) == 0 ||  # special case of empty ComponentVector
        return return (; isvalid=false,
            msg=String127("Expected paropt.par <: ComponentVector, but was not."))
    u0c = attach_axis((1:axis_length(ax_state_scalar)), ax_state_scalar)
    for k in keys(paropt.state)
        k ∉ keys(u0c) && return (; isvalid=false,
            msg=String127("Expected optimized paropt.state.$k to be part of state, " *
                          "but was not."))
        length(paropt.state[k]) != length(u0c[k]) && return (; isvalid=false,
            msg=String127("Expected optimized paropt.state.$k to be of length " *
                          "$(length(u0c[k])) but had length $(length(paropt.state[k]))"))
    end
    pc = attach_axis((1:axis_length(ax_par)), ax_par)
    for k in keys(paropt.par)
        length_par_k = length(paropt.par[k])::Int
        k ∉ keys(pc) && return (; isvalid=false,
            msg=String127("Expected optimized paropt.par.$k to be part of parameters, " *
                          "but was not."))
        length_par_k != length(pc[k]) && return (; isvalid=false,
            msg=String127("Expected optimized paropt.par.$k to be of length " *
                          "$(length(pc[k])) but had length $length_par_k"))
    end
    return (; isvalid=true, msg=String127(""))
end

"""
    vcat_statesfirst(cvs...; system)

Concatenate ComponentVectors, `cvs`, but move state entries before parameter entries.   
"""
function vcat_statesfirst(cvs...; system::AbstractSystem)
    popt0 = reduce(vcat, cvs)
    pset = LoggingExtras.withlevel(Logging.Error) do
        pset = ODEProblemParSetter(system, popt0)
    end
    popt_template = label_paropt_flat1(pset, 1:count_paropt(pset))
    popt0[keys(popt_template)]
end

