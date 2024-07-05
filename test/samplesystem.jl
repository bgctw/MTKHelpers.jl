
using ModelingToolkit: t_nounits as t, D_nounits as D
using DomainSets
using Memoize

"""
    samplesystem(; name, τ = 3.0, i = 0.1, p1 = 1.1, p2 = 1.2)

Defines a simple ODESystem of exponential decay to p1/p2 with rate τ
"""
@memoize function samplesystem(; name, τ = 3.0, i = 0.1, p1 = 1.1, p2 = 1.2)
    sts = Symbolics.@variables x(t) RHS(t)  # RHS is observed
    ps = @parameters τ=τ p1=p1 p2=p2 i=i       # parameters
    ODESystem([
            RHS ~ i - p1 * x^2 + (p2 - x) / τ,
            D(x) ~ RHS
        ], t, sts, ps; name)
end

@memoize function samplesystem_const(RHS0; name)
    # demonstrating override_system by setting the RHS to constant first order rate
    m = samplesystem(; name)
    @unpack RHS, x = m
    ps = @parameters RHS_0 = RHS0
    eqs = [RHS ~ RHS_0 * x]
    sys_ext = override_system(eqs, m; name, ps)
end

# function samplesystem_vec(; name, τ = 3.0, i=0.1, p = [1.1, 1.2])
#     ncomp = 2
#     sts = @variables x(t)[1:2] dx(t)[1:2]  # observed dx now can be accessed
#     ps = @parameters τ=τ p p2=p[2] i=i       # parameters
#     ODESystem([
#             dx₁(t) ~ i - p1 * x₁(t)^2 + (p2 - x₁(t)) / τ,
#             dx₂(t) ~ i - p1 * x₂(t)^2 + (p2 - x₂(t)) / τ,
#             D(x₁) ~ dx₁(t),
#             D(x₂) ~ dx₂(t),
#         ], t, sts, ps; name)
# end

@memoize function get_sys_ex_scalar()
    # states and parameters are single entries
    u1 = CA.ComponentVector(L = 10.0)
    p1 = CA.ComponentVector(k_L = 1.0, k_R = 1 / 20, m = 2.0)
    popt1 = CA.ComponentVector(L = 10.1, k_L = 1.1, k_R = 1 / 20.1)
    popt1s = CA.ComponentVector(state = (L = 10.1,), par = (k_L = 1.1, k_R = 1 / 20.1))
    function get_sys1()
        sts = @variables L(t)
        ps = @parameters k_L, k_R, m
        eq = [D(L) ~ 0]
        sys1 = ODESystem(eq, t, sts, vcat(ps...); name = :sys1)
    end
    sys1 = structural_simplify(get_sys1())
    prob_sys1 = ODEProblem(
        sys1, get_system_symbol_dict(sys1, u1), (0.0, 1.0),
        get_system_symbol_dict(sys1, p1))
    (; u0 = u1, p = p1, popt = popt1s, prob = prob_sys1)
end

@memoize function get_sys_ex_vec()
    # entries with substructure
    u1c = CA.ComponentVector(a = [1, 21, 22.0])
    #p1c = CA.ComponentVector(b = (b1 = 0.1, b2 = 0.2), c = [0.01, 0.02], d = 3.0)
    p1c = CA.ComponentVector(b = [0.1, 0.2], c = [0.01, 0.02], d = 3.0)
    # as long as ComponentArrays does not support Axis-indexing, focus on top-level components rather than implementing this indexing in MTKHelpers
    # note: no a1 no b, 
    # a2 and c need to have correct length for updating
    #poptc = CA.ComponentVector(a=(a2=1:2,), c=1:2) 
    poptc = vcat(u1c[CA.KeepIndex(:a)], p1c[(:b, :c)])
    poptcs = CA.ComponentVector(state = u1c[CA.KeepIndex(:a)], par = p1c[(:b, :c)])
    function get_sys2()
        @variables a(..)[1:3]
        sts = vcat([a(t)[i] for i in 1:3])
        ps = @parameters b[1:2], c[1:2], d
        eq = vcat([D(a(t)[i]) ~ 0 for i in 1:3])
        sys2 = ODESystem(eq, t, sts, vcat(ps...); name = :sys2)
    end
    sys2 = structural_simplify(get_sys2())
    prob_sys2 = ODEProblem(
        sys2, get_system_symbol_dict(sys2, u1c), (0.0, 1.0),
        get_system_symbol_dict(sys2, p1c))
    (; u0 = u1c, p = p1c, popt = poptcs, prob = prob_sys2)
end

@memoize function get_sys_ex_pde()
    @parameters t z
    z_m = +0.3 # maximum depth in m, change to positive depth
    n_z = 16
    #z_grid = collect(range(0,z_m, length=n_z))
    z_grid = grid_exp(n_z, z_m, 3.0)
    dz = z_m / (n_z - 1) # control of same grid represented by a single number
    discretization = MOLFiniteDifference([z => z_grid], t;
        advection_scheme = UpwindScheme(), approx_order = 2)

    @parameters k_Y Y0 i_Y ω i_Y_agr[1:2]
    @variables Y(..) i_Yo(..) adv_Yo(..) dec_Y(..) Y_t(..) Y_zz(..) Yi(..) i_Yi(..) Y_z(..) adv_Yi(..)
    ∂_t = Differential(t)
    ∂_z = Differential(z)

    z_min = 0.0  # directly using 0.0 in Integral causes error in solution wrapping
    Iz = Integral(z in DomainSets.ClosedInterval(z_min, z))

    eqs0 = [
        ∂_t(Y(t, z)) ~ Y_t(t, z),
        Y_t(t, z) ~ i_Yo(t, z) - dec_Y(t, z) + adv_Yo(t, z),
        #i_Yo(t, z) ~ i_Yz(t, z, i_Y),         # observable input Y specified below
        dec_Y(t, z) ~ k_Y * Y(t, z),          # observable of decomposition 
        #adv_Yo(t, z) ~ -ω * ∂_z(Y(t, z)),  # observable advective flux of Y
        adv_Yo(t, z) ~ -ω * ∂_z(Y(t, z)),  # observable advective flux of Y
        # further observables that are not used in eqs
        Y_z(t, z) ~ ∂_z(Y(t, z)),
        Yi(t, z) ~ Iz(Y(t, z)),
        i_Yi(t, z) ~ Iz(i_Yo(t, z)), # integral over litter inputs to check
        adv_Yi(t, z) ~ Iz(adv_Yo(t, z)) # integral over advection inputs
    ]

    i_Yz(t, z, i_Y) = Dz_exp(z, z_m, 4.0) * i_Y
    #@register_symbolic i_Yz(t, z, i_Y)

    eqs = vcat(eqs0, [
        i_Yo(t, z) ~ i_Yz(t, z, i_Y)         # observable input Y
    ])

    fgauss(x, μ, σ2) = 1 / sqrt(2 * pi * σ2) * exp(-(x - μ)^2 / (2 * σ2))
    fagr(t, i_Y_agr, i_Y_agr_pulse) = i_Y_agr + i_Y_agr_pulse * fgauss(t, 80, 2^2)
    #@register_symbolic fagr(t, i_Y_agr, i_Y_agr_pulse)

    bcs = [
        Y(0, z) ~ Dz_lin(z, z_m) * Y0, # initial constant with depth,
        ω * Y(t, 0) ~ fagr(t, i_Y_agr[1], i_Y_agr[2]), # specified flux at upper boundary
        ∂_z(Y(t, z_m)) ~ 0, # negligible change in concentration at lower boundary
        #
        # following are only necessary because observables boundaries 
        # need to be specified with vector grid
        ∂_z(dec_Y(t, 0)) ~ 0.0,
        ∂_z(dec_Y(t, z_m)) ~ 0.0,
        ∂_z(i_Yo(t, 0)) ~ 0.0,
        ∂_z(i_Yo(t, z_m)) ~ 0.0,
        #∂_z(adv_Yo(t, 0)) ~ 0.0, 
        ω * adv_Yo(t, 0) ~ fagr(t, i_Y_agr[1], i_Y_agr[2]),
        ∂_z(adv_Yo(t, z_m)) ~ 0.0,
        ∂_z(Y_t(t, 0)) ~ 0.0,
        ∂_z(Y_t(t, z_m)) ~ 0.0,
        ∂_z(Yi(t, 0)) ~ 0.0,
        ∂_z(Yi(t, z_m)) ~ 0.0,
        ∂_z(i_Yi(t, 0)) ~ 0.0,
        ∂_z(i_Yi(t, z_m)) ~ 0.0,
        ∂_z(Y_z(t, 0)) ~ 0.0,
        ∂_z(Y_z(t, z_m)) ~ 0.0,
        ∂_z(adv_Yi(t, 0)) ~ 0.0,
        ∂_z(adv_Yi(t, z_m)) ~ 0.0
    ]

    # Space and time domains
    domains = [t ∈ Interval(0.0, 500.0), z ∈ Interval(0.0, z_m)]

    # PDE system
    state_vars = [Y(t, z), i_Yo(t, z), adv_Yo(t, z), dec_Y(t, z),
        Y_t(t, z), Yi(t, z), i_Yi(t, z), Y_z(t, z), adv_Yi(t, z)]
    params = [
        k_Y => 2.0,
        Y0 => 200.0,
        i_Y => 50.0,
        ω => 0.01,
        i_Y_agr[1] => 10.0, # base input
        i_Y_agr[2] => 50.0 # pulse at t=80       
    ]
    @named pdesys = PDESystem(eqs, bcs, domains, [t, z], state_vars, params)
    prob = discretize(pdesys, discretization)

    p = CA.ComponentVector(k_Y = 1/10, Y0 = 80.0, i_Y = 10.0,  i_Y_agr = [2.0, 50.0]) 
    par_new = p[[:Y0, :i_Y, :i_Y_agr]] 
    #
    zs = get_1d_grid(get_system(prob))
    state_pos = get_1d_state_pos(get_system(prob))
    u0 = CA.ComponentVector(Y = Dz_lin.(zs[state_pos], z_m) * par_new.Y0)
    
    paropt = CA.ComponentVector(state = u0, par = par_new)
    (; u0, p = par_new, popt = paropt, prob, zs, state_pos, pdesys)
end