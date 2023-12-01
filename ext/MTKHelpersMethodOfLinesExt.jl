module MTKHelpersMethodOfLinesExt

function __init__()
    @info "MTKHelpers: loading MTKHelpersMethodOfLinesExt"
end

isdefined(Base, :get_extension) ? (using MethodOfLines) : (using ..MethodOfLines)
using MTKHelpers

using ModelingToolkit
using ModelingToolkit: AbstractODESystem, AbstractSystem
using SciMLBase: SciMLBase, AbstractODEProblem
using Chain
using DomainSets # for the example system
using ComponentArrays

function MTKHelpers.get_discrete_space(prob::AbstractODEProblem)
    sys = get_system(prob)
    get_discrete_space(sys)
end
function MTKHelpers.get_discrete_space(sys::AbstractSystem)
    hasfield(typeof(ModelingToolkit.get_metadata(sys)), :discretespace) ||
        error("Cannot get ConcreteSpace for the system associated to this problem")
    ModelingToolkit.get_metadata(sys).discretespace
end

function MTKHelpers.get_1d_grid(prob::AbstractODEProblem)
    ds = get_discrete_space(prob)
    length(ds.axies) == 1 ||
        error("Cannot get 1d grid of a problem is associated with a system of " *
              string(length(ds.axies)) * " spatial variables.")
    first(ds.axies).second
end

function MTKHelpers.get_1d_state_pos(prob::AbstractODEProblem)
    sys = MTKHelpers.get_system(prob)
    ds = MTKHelpers.get_discrete_space(sys)
    length(ds.axies) == 1 ||
        error("Cannot get 1d grid of a problem is associated with a system of " *
              string(length(ds.axies)) * " spatial variables.")
    # ds.discvars relates to non-simplified and holds all observables -> need sys
    us = MTKHelpers.get_states(sys)
    # get_gridloc  does not work for Symbolic-Array-States
    #g_vars = MethodOfLines.get_gridloc.(us, (ds,)) # .Tuple{Num, d_vector_pos}
    # @chain g_vars begin
    #     # only those discretized states that are the same as the first (same pde state)
    #     filter(x -> Symbol(first(x)) == Symbol(first(g_vars[1])), _)
    #     getindex.(2)  # extract the position vector
    #     first.()      # take the first dimension-entry
    # end
    z_grid = get_1d_grid(prob)
    @chain us begin
        # only getindex vars, i.e. discretized ones
        filter(x -> isequal(operation(x), getindex), _)
        # only the first one of these (having the same 1st argument - the name)
        filter(x -> isequal(arguments(x)[1], arguments(first(_))[1]), _)
        # extract the position (second argument)
        getindex.(arguments.(_), 2)
        # index into z_grid
        # getindex.(Ref(z_grid), _)
    end
end

function MTKHelpers.example_pde_problem(; name = "mp")
    @parameters t z
    z_m = +0.3 # maximum depth in m, change to positive depth
    n_z = 16
    z_grid = grid_exp(n_z, z_m, 3) # grid denser at low depth (top)
    #z_grid = collect(range(0,z_m, length=n_z))
    discretization_grid = MOLFiniteDifference([z => z_grid], t;
        advection_scheme = UpwindScheme(), approx_order = 2)
    # dzs = diff(z_grid)
    # dzsl = vcat(dzs[1]/2, dzs[2:end], dzs[end]/2) # assume first and last layer only half 
    @parameters k_Y Y0 i_Y ω i_Y_agr[1:2]
    @variables Y(..) i_Yo(..) adv_Yo(..) dec_Y(..) Y_t(..) Y_zz(..) Yi(..) i_Yi(..) Y_z(..) adv_Yi(..)
    ∂_t = Differential(t)
    ∂_z = Differential(z)
    params = [
        k_Y => 2.0,
        Y0 => 200.0,
        i_Y => 50.0,
        ω => 0.01,
        i_Y_agr[1] => 10.0, # base input
        i_Y_agr[2] => 50.0, # pulse at t=80       
    ]
    #    
    # below-ground inputs (density across depth) exponentially decreasing with depth
    i_Yz(t, z, i_Y) = Dz_exp(z, z_m, 4.0) * i_Y  # inputs that integrate to i_Y independent of time
    #@register_symbolic i_Yz(t, z, i_Y)
    #
    # above-ground inputs constant with a pulse at t=80
    fgauss(x, μ, σ2) = 1 / sqrt(2 * pi * σ2) * exp(-(x - μ)^2 / (2 * σ2))
    fagr(t, i_Y_agr, i_Y_agr_pulse) = i_Y_agr + i_Y_agr_pulse * fgauss(t, 80, 2^2)
    # @register_symbolic fagr(t, i_Y_agr, i_Y_agr_pulse)
    #
    z_min = 0.0  # directly using 0.0 in Integral causes error in solution wrapping
    Iz = Integral(z in DomainSets.ClosedInterval(z_min, z))
    #
    # Space and time domains
    domains = [t ∈ Interval(0.0, 500.0), z ∈ Interval(0.0, z_m)]
    #
    eqs2 = [
        ∂_t(Y(t, z)) ~ i_Yz(t, z, i_Y) - dec_Y(t, z) + adv_Yo(t, z),
        dec_Y(t, z) ~ k_Y * Y(t, z),          # observable of decomposition 
        adv_Yo(t, z) ~ -ω * ∂_z(Y(t, z)),     # observable advective flux of Y
        Yi(t, z) ~ Iz(Y(t, z)),               # observable integral from 0 to depth z
    ]
    bcs2 = [
        Y(0, z) ~ 1 / z_m * Y0,
        ω * Y(t, 0) ~ fagr(t, i_Y_agr[1], i_Y_agr[2]),
        # only needed for observables being discretized
        ∂_z(Y(t, z_m)) ~ 0.0,
        adv_Yo(t, 0) ~ fagr(t, i_Y_agr[1], i_Y_agr[2]),
        ∂_z(adv_Yo(t, z_m)) ~ 0.0,
        ∂_z(dec_Y(t, 0)) ~ 0.0,
        ∂_z(dec_Y(t, z_m)) ~ 0.0,
        ∂_z(Yi(t, 0)) ~ 0.0,
        ∂_z(Yi(t, z_m)) ~ 0.0,
    ]
    state_vars2 = [Y(t, z), dec_Y(t, z), adv_Yo(t, z), Yi(t, z)]
    @named pdesys2 = PDESystem(eqs2, bcs2, domains, [t, z], state_vars2, params)
    # Convert the PDE problem into an ODE problem
    prob2 = discretize(pdesys2, discretization_grid)
end

function MTKHelpers.example_pde_problem_arrstate(; name = "ma")
    @parameters t z
    z_m = +0.3 # maximum depth in m, change to positive depth
    n_z = 16
    z_grid = collect(range(0, z_m, length = n_z))
    dz = z_m / (n_z - 1) # control of same grid represented by a single number
    discretization = MOLFiniteDifference([z => z_grid], t;
        advection_scheme = UpwindScheme(), approx_order = 2)
    #    
    n_comp = 2
    @parameters k_Y Y0 ω
    @variables Y(..)[1:n_comp] adv_Yo(..)[1:n_comp] dec_Y(..)[1:n_comp] Y_t(..)[1:n_comp]
    @variables Yi(..)[1:n_comp]
    ∂_t = Differential(t)
    ∂_z = Differential(z)
    params = [
        k_Y => 2.0,
        Y0 => 200.0,
        ω => 0.01,
    ]
    #
    z_min = 0.0  # directly using 0.0 in Integral causes error in solution wrapping
    Iz = Integral(z in DomainSets.ClosedInterval(z_min, z))
    #
    eqs_s = [
        Y_t(t, z)[1] ~ -dec_Y(t, z)[1] + adv_Yo(t, z)[1] + 20,
        Y_t(t, z)[2] ~ -dec_Y(t, z)[2] + adv_Yo(t, z)[2],
    ]
    eqs_i = [[
        ∂_t(Y(t, z)[i]) ~ Y_t(t, z)[i],
        dec_Y(t, z)[i] ~ k_Y * Y(t, z)[i],          # observable of decomposition 
        adv_Yo(t, z)[i] ~ -ω * ∂_z(Y(t, z)[i]),  # observable advective flux of Y
        # further observables that are not used in eqs
        Yi(t, z)[i] ~ Iz(Y(t, z)[i]),
    ] for i in 1:n_comp]
    eqs = vcat(eqs_s, eqs_i...)
    #
    bcs_s = [
    ]
    bcs_i = [[
        #Y(0, z) ~ Dz_exp(z, z_m, 2.0) * Y0, # initial exponential distribution with depth
        Y(0, z)[i] ~ Y0 / z_m, #Dz_lin(z, z_m) * Y0, # initial constant with depth, modified by set problem.u0
        ω * Y(t, 0)[i] ~ 10.0, # specified flux at upper boundary
        ∂_z(Y(t, z_m)[i]) ~ 0, # negligible change in concentration at lower boundary
        #following are only necessary because observables need to be specified with vector grid
        ∂_z(dec_Y(t, 0)[i]) ~ 0.0,
        ∂_z(dec_Y(t, z_m)[i]) ~ 0.0,
        #∂_z(adv_Yo(t, 0)[i]) ~ 0.0, 
        ω * adv_Yo(t, 0)[i] ~ 10.0,
        ∂_z(adv_Yo(t, z_m)[i]) ~ 0.0,
        ∂_z(Y_t(t, 0)[i]) ~ 0.0,
        ∂_z(Y_t(t, z_m)[i]) ~ 0.0,
        ∂_z(Yi(t, 0)[i]) ~ 0.0,
        ∂_z(Yi(t, z_m)[i]) ~ 0.0,
    ] for i in 1:n_comp]
    bcs = vcat(bcs_s, bcs_i...)
    #
    # Space and time domains
    domains = [t ∈ Interval(0.0, 500.0), z ∈ Interval(0.0, z_m)]
    #
    # PDE system
    state_vars = [Y(t, z), adv_Yo(t, z), dec_Y(t, z), Y_t(t, z), Yi(t, z)]
    @named pdesys = PDESystem(eqs, bcs, domains, [t, z], state_vars, params)
    prob = discretize(pdesys, discretization)
end

end
