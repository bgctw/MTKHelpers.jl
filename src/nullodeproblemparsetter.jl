# simplified version of ODEProblemParSetter with empty paropt
# specialized for increasing performance
# does not modify Problem on remake

struct NullODEProblemParSetter <: AbstractODEProblemParSetter
    ax_state::AbstractAxis
    ax_state_scalar::AbstractAxis
    ax_par::AbstractAxis
    par_ind::AbstractVector # propb.ps -> vector
end

struct NullODEProblemParSetterConcrete{
    SA <: AbstractAxis, PA <: AbstractAxis, SSA <: AbstractAxis, NPS, EIP} <: AbstractODEProblemParSetter
    ax_state::SA
    ax_state_scalar::SSA
    ax_par::PA
    par_ind::SVector{NPS, EIP} # propb.ps -> vector
end

function get_concrete(pset::NullODEProblemParSetter)
    NullODEProblemParSetterConcrete(
        pset.ax_state, pset.ax_state_scalar, pset.ax_par, SVector{length(pset.par_ind)}(pset.par_ind)
    )
end


function NullODEProblemParSetter(state_template, par_template, system::AbstractODESystem)
    ax_state = _get_axis(state_template)
    ax_par = _get_axis(par_template)
    scalar_num_map = get_scalar_num_map(system)
    (ax_state_scalar, ax_paropt_scalar) = scalarize_par_and_paroptstate(ax_state, CA.Axis(state = 1:0, par = 1:0), scalar_num_map)
    _dict_nums = get_system_symbol_dict(system)
    par_nums = axis_length(ax_par) == 0 ? eltype(values(_dict_nums))[] : vcat(
        Symbolics.scalarize.(
            getindex.(Ref(_dict_nums), keys(ax_par)))...)
    par_ind = SII.parameter_index.(Ref(system), par_nums)
    NullODEProblemParSetter(ax_state, ax_state_scalar, ax_par, par_ind)
end

function NullODEProblemParSetter(sys::AbstractODESystem)
    NullODEProblemParSetter(
        #Axis(symbol_op_scalar.(unknowns(sys))),
        axis_of_nums(unknowns(sys)),
        axis_of_nums(parameters(sys)),
        sys)
end


NullODEProblemParSetterU = Union{NullODEProblemParSetter, NullODEProblemParSetterConcrete}

axis_state(ps::NullODEProblemParSetterU) = ps.ax_state
axis_state_scalar(ps::NullODEProblemParSetterU) = ps.ax_state_scalar
axis_par(ps::NullODEProblemParSetterU) = ps.ax_par
axis_paropt(ps::NullODEProblemParSetterU) = CA.Axis(state = 1:0, par = 1:0)
axis_paropt_scalar(ps::NullODEProblemParSetterU) = CA.Axis(state = 1:0, par = 1:0)
axis_paropt_flat1(ps::NullODEProblemParSetterU) = FlatAxis()

classes_paropt(::NullODEProblemParSetterU) = (:state, :par)

function get_paropt(pset::NullODEProblemParSetterU, prob::AbstractODEProblem; kwargs...)
    typeof(promote())
    T = promote_type(eltype(prob.u0), eltype(prob.ps[pset.par_ind]))
    T[]
end

function get_par(pset::NullODEProblemParSetterU, prob::SciMLBase.AbstractSciMLProblem)
    prob.ps[pset.par_ind]
end

# override warning on non-flat version
function label_paropt_flat1(pset::NullODEProblemParSetterU, popt; omitwarning = false)
    ComponentVector{eltype(popt)}()
end

# fast version: Problem not modified
function remake_pset(prob::AbstractODEProblem, popt, pset::NullODEProblemParSetterU)
    prob
end



