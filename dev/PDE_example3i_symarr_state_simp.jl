using OrdinaryDiffEq, ModelingToolkit, DomainSets
using MethodOfLines
#using MTKHelpers

@parameters t z
z_m = +0.3 # maximum depth in m, change to positive depth
n_z = 16
z_grid = collect(range(0, z_m, length = n_z))
dz = z_m / (n_z - 1) # control of same grid represented by a single number
discretization = MOLFiniteDifference([z => z_grid], t;
    advection_scheme = UpwindScheme(), approx_order = 2)

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

z_min = 0.0  # directly using 0.0 in Integral causes error in solution wrapping
Iz = Integral(z in DomainSets.ClosedInterval(z_min, z))

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

# Space and time domains
domains = [t ∈ Interval(0.0, 500.0), z ∈ Interval(0.0, z_m)]

# PDE system
state_vars = [Y(t, z), adv_Yo(t, z), dec_Y(t, z), Y_t(t, z), Yi(t, z)]
@named pdesys = PDESystem(eqs, bcs, domains, [t, z], state_vars, params)
# Convert the PDE problem into an ODE problem
prob = discretize(pdesys, discretization)
