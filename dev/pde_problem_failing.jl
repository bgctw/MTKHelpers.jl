Pkg.activate(; temp=true)
Pkg.add(["MethodOfLines","ModelingToolkit","DomainSets"]) 

using MethodOfLines, ModelingToolkit
using DomainSets
#using MTKHelpers

"""
    Dz_exp(x,x_m,b)
    Iz_exp(x,x_m,b)

Density function decreasing exponentially from zero with e-folding time b 
``x \\in [0,x_m]`` for ``x_m > 0`` for which the intral ``\\int_0^{x_m}  = 1``.
"""
function Dz_exp(x, x_m, b)
    b / (1 - exp(-b * x_m)) * exp(-b * x)
end,
function Iz_exp(x, x_m, b)
    1 / (exp(-b * x_m) - 1) * (exp(-b * x) - 1)
end


    @parameters t z
    z_m = +0.3 # maximum depth in m, change to positive depth
    n_z = 4#16
    #z_grid = grid_exp(n_z, z_m, 3) # grid denser at low depth (top)
    z_grid = [0.0, 0.003480165515666162, 0.007730849275354716, 0.012922646143504952, 0.019263921158071014, 0.027009171951114137, 0.036469242632379266, 0.04802379905486663, 0.06213656613861004, 0.07937393877996625, 0.10042771326755343, 0.1261428514963746, 0.1575513922555268, 0.19591387056854145, 0.24276990738991716, 0.3]
    #z_grid = collect(range(0,z_m, length=n_z))
    discretization_grid = MOLFiniteDifference([z => z_grid], t;
        advection_scheme = UpwindScheme(), approx_order = 2)
    # dzs = diff(z_grid)
    # dzsl = vcat(dzs[1]/2, dzs[2:end], dzs[end]/2) # assume first and last layer only half 
    @parameters k_Y Y0 i_Y ω i_Y_agr[1:2] X0
    @variables Y(..) i_Yo(..) adv_Yo(..) dec_Y(..) Y_t(..) Y_zz(..) Yi(..) i_Yi(..) Y_z(..)
    @variables adv_Yi(..) X(..)
    ∂_t = Differential(t)
    ∂_z = Differential(z)
    params = [
        k_Y => 2.0,
        Y0 => 200.0,
        i_Y => 50.0,
        ω => 0.01,
        i_Y_agr[1] => 10.0, # base input
        i_Y_agr[2] => 50.0, # pulse at t=80   
        X0 => 2.0,
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
        ∂_t(X(t)) ~ -k_Y * X(t),              # state variable that is not discretized
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
        X(0) ~ X0,
    ]
    state_vars2 = [Y(t, z), dec_Y(t, z), adv_Yo(t, z), Yi(t, z), X(t)]
    @named pdesys2 = PDESystem(eqs2, bcs2, domains, [t, z], state_vars2, params)
    # Convert the PDE problem into an ODE problem
    prob2 = discretize(pdesys2, discretization_grid)
