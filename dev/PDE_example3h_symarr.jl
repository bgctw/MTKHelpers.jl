using OrdinaryDiffEq, ModelingToolkit, DomainSets
using MethodOfLines
using MTKHelpers

@parameters t z
z_m = +0.3 # maximum depth in m, change to positive depth
n_z = 16
#z_grid = collect(range(0,z_m, length=n_z))
z_grid = grid_exp(n_z, z_m, 3.0)
#dz = z_m/(n_z-1) # control of same grid represented by a single number
#z_grid = z_m/(n_z-1)
#scatter(z_grid, yflip=true)
discretization = MOLFiniteDifference([z => z_grid], t;
    advection_scheme = UpwindScheme(), approx_order = 2)
#advection_scheme=WENOScheme(), approx_order = 2)
dzs = diff(z_grid)
dzsl = vcat(dzs[1] / 2, dzs[2:end], dzs[end] / 2) # assume first and last layer only half 

@parameters k_Y Y0 i_Y ω i_Y_agr[1:2]
@variables Y(..) i_Yo(..) adv_Yo(..) dec_Y(..) Y_t(..) Y_zz(..) Yi(..) i_Yi(..) Y_z(..) adv_Yi(..)
∂_t = Differential(t)
∂_z = Differential(z)
∂_zz = Differential(z)^2
∂_zzz = Differential(z)^3
params = [
    k_Y => 2.0,
    Y0 => 200.0,
    i_Y => 50.0,
    #ω => -0.01,
    ω => 0.01,
    i_Y_agr[1] => 10.0, # base input
    i_Y_agr[2] => 50.0, # pulse at t=80       
]

#i_Yz(t,z,i_Y) = Dz_lin(z) * i_Y   # inputs that integrate to i_Y independent of time
i_Yz(t, z, i_Y) = Dz_exp(z, z_m, 4.0) * i_Y  # inputs that integrate to i_Y independent of time
#@register_symbolic i_Yz(t, z, i_Y)

tmp_f = () -> begin
    _iY = i_Yz.(0, z_grid, 40.0)
    sum(_iY .* dzsl) # coarse approximation of an integral
    scatter(_iY, z_grid, yflip = true, xlabel = "iY (g/m/s)")
end

fgauss(x, μ, σ2) = 1 / sqrt(2 * pi * σ2) * exp(-(x - μ)^2 / (2 * σ2))
fagr(t, i_Y_agr, i_Y_agr_pulse) = i_Y_agr + i_Y_agr_pulse * fgauss(t, 80, 2^2)
#@register_symbolic fagr(t, i_Y_agr, i_Y_agr_pulse)

tmpf = () -> begin
    x = 0:200
    plot(x, 20 * fgauss.(x, 80, 10^2))
    plot(x, fagr.(x, 2.0, 2.0 * 10))
    dx = 0.02
    isapprox(sum(fagr.(70:dx:90, 0, 50) * dx), 50, atol = 1e-3) # integral of agr_pulse matches
end

z_min = 0.0  # directly using 0.0 in Integral causes error in solution wrapping
Iz = Integral(z in DomainSets.ClosedInterval(z_min, z))

eqs = [
    ∂_t(Y(t, z)) ~ Y_t(t, z),
    Y_t(t, z) ~ i_Yo(t, z) - dec_Y(t, z) + adv_Yo(t, z),
    i_Yo(t, z) ~ i_Yz(t, z, i_Y),         # observable input Y
    dec_Y(t, z) ~ k_Y * Y(t, z),          # observable of decomposition 
    #adv_Yo(t, z) ~ -ω * ∂_z(Y(t, z)),  # observable advective flux of Y
    adv_Yo(t, z) ~ -ω * ∂_z(Y(t, z)),  # observable advective flux of Y
    # further observables that are not used in eqs
    Y_z(t, z) ~ ∂_z(Y(t, z)),
    #Y_zz(t,z) ~ ∂_zz(Y(t, z)),
    Yi(t, z) ~ Iz(Y(t, z)),
    i_Yi(t, z) ~ Iz(i_Yo(t, z)), # integral over litter inputs to check
    adv_Yi(t, z) ~ Iz(adv_Yo(t, z)), # integral over advection inputs
]

bcs = [
    #Y(0, z) ~ Dz_exp(z, z_m, 2.0) * Y0, # initial exponential distribution with depth
    Y(0, z) ~ Dz_lin(z, z_m) * Y0, # initial constant with depth, modified by set problem.u0
    ω * Y(t, 0) ~ fagr(t, i_Y_agr[1], i_Y_agr[2]), # specified flux at upper boundary
    ∂_z(Y(t, z_m)) ~ 0, # negligible change in concentration at lower boundary
    #
    # following are only necessary because observables need to be specified with vector grid
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
    ∂_z(adv_Yi(t, z_m)) ~ 0.0,
]

# Space and time domains
domains = [t ∈ Interval(0.0, 500.0), z ∈ Interval(0.0, z_m)]

# PDE system
state_vars = [
    Y(t, z),
    i_Yo(t, z),
    adv_Yo(t, z),
    dec_Y(t, z),
    Y_t(t, z),
    Yi(t, z),
    i_Yi(t, z),
    Y_z(t, z),
    adv_Yi(t, z),
]
#state_vars = [Y(t, z), i_Yo(t,z), adv_Yo(t,z), dec_Y(t,z), Y_t(t,z), Y_z(t,z),]
#state_vars = [Y(t, z), i_Yo(t,z), adv_Yo(t,z), dec_Y(t,z), Y_t(t,z), Y_z(t,z), Yi(t,z),]
#state_vars = [Y(t, z), i_Yo(t,z), adv_Yo(t,z), dec_Y(t,z), Y_t(t,z)]
@named pdesys = PDESystem(eqs, bcs, domains, [t, z], state_vars, params)
# Convert the PDE problem into an ODE problem
prob = discretize(pdesys, discretization)

tmp_f = () -> begin
    # get the ODESystem to construct ProbParameterSetter
    odesys_full, tspan = symbolic_discretize(pdesys, discretization)
    #odesys_full, tspan = symbolic_discretize(pdesys2, discretization);
    odesys = structural_simplify(odesys_full)
    observed(odesys)
    equations(odesys)
    parameters(odesys)
    unknowns(odesys)
    # removes observables - state only contains Y
end

# Solve ODE problem
#sol = solve(prob, Tsit5(), saveat=0.2);
#sol = solve(prob, Tsit5());
sol = solve(prob, TRBDF2());
#sol = solve(prob2, TRBDF2());
#plot(sol[dec_Y(t,z)][1,:])
#plot(sol[dec_Y(t,z)][end,:])
#plot(sol[Y(t,z)][1,:]); plot!(sol[Y(t,z)][end,:])
#plot(sol[adv_Yo(t,z)][1,:]); plot!(sol[adv_Yo(t,z)][end,:])

po = prob.p;
po[4] = 0.2; # higher leaching/advection rate
probo = remake(prob, p = po)
solp = solo = solve(probo, TRBDF2());

# compare end proviles between diffusion rates
# z_m:dz:0 (top is displayed last)
hcat(sol[Y(t, z)][end, :], solo[Y(t, z)][end, :]) # z_m:dz:0 (top is displayed last)

# last time derivates
solp[t][(end .- 5):end, :]
solp[Y_t(t, z)][(end .- 5):end, :]
# should be near zero at the end, except for top -balanced by agr input
solp[adv_Yo(t, z)][(end .- 5):end, :]
solp[Y(t, z)][(end .- 5):end, :] # in steady state

# advective fluxes are of magnitude of i_Yagr input -> same magnitude as i_Y
# different at depth
hcat(solp[i_Yo(t, z)][end, :], -solp[dec_Y(t, z)][end, :], solp[adv_Yo(t, z)][end, :])
# input == -dec + adv    , except top layer (last row where there is additional agr input)
hcat(solp[i_Yo(t, z)][end, :], -solp[dec_Y(t, z)][end, :] + solp[adv_Yo(t, z)][end, :])

using Plots
#plot(sol[Y(t,z)][1,:], discrete_z, label="ω=0.01");
plot(sol[Y(t, z)][end, :], discrete_z, label = "ω=0.01");
plot!(solo[Y(t, z)][end, :], discrete_z, label = "ω=0.2")

#solp = sol;
#solp = solo;
plot(solp[i_Yo(t, z)][end, :] ./ maximum(solp[i_Yo(t, z)][end, :]),
    discrete_z,
    label = "∂i_z",
    xlims = (0, +Inf));
plot!(solp[Y(t, z)][end, :] ./ maximum(solp[Y(t, z)][end, :]), discrete_z, label = "∂Y_z")
#plot!(solp[dec_Y(t,z)][end,:] ./ maximum(solp[dec_Y(t,z)][end,:]), discrete_z, label="∂decY_z");
plot!(solp[adv_Yo(t, z)][end, :] ./ maximum(solp[adv_Yo(t, z)][end, :]),
    discrete_z,
    label = "∂advY_z")

# i + adv balances dec, except at top layer where there is additional advective influx
plot(solp[i_Yo(t, z)][end, :] .* dz, discrete_z, label = "∂i_z", xlims = (0, +Inf));
plot!(solp[dec_Y(t, z)][end, :] .* dz, discrete_z, label = "∂decY_z");
plot!(solp[adv_Yo(t, z)][end, :] .* dz, discrete_z, label = "∂advY_z")

# First layer
layer = 1
plot(solp[t], solp[Y(t, z)][:, end + 1 - layer])
plot(solp[t], solp[adv_Yo(t, z)][:, end + 1 - layer])

# integrated litter input matches parameter? 
isapprox(sol[i_Yi(t, z)][end, end], Dict(params)[i_Y]; rtol = 0.01)
isapprox(solp[i_Yi(t, z)][end, end], par_new.i_Y; rtol = 0.01)

# Yi (integrated across profile) vs. time
plot(discrete_to,
    solo[Yi(t, z)][:, end],
    label = "ω=0.01",
    xlim = (0, 10),
    xlab = "Time (yr)",
    ylab = "Y (g)");
plot!(discrete_t, sol[Yi(t, z)][:, end], label = "ω=0.2")

solo[∂_z(Y(t, z))]

#---------------- pulse experiment -----------------------------
using ComponentArrays
# only advection: low decomposition and low below ground input
# without any abr inputs
par_new = ComponentVector(k_Y = 1e-6, Y0 = 100.0, i_Y = 1e-12,
    ω = 0.01, i_Y_agr = [0.0, 0.0])
# with abr input rate
par_new = ComponentVector(k_Y = 1e-6, Y0 = 100.0, i_Y = 1e-12,
    ω = 0.01, i_Y_agr = [10.0, 0.0])
# with a pulse no input rate
par_new = ComponentVector(k_Y = 1e-6, Y0 = 0.0, i_Y = 1e-12,
    ω = 0.01, i_Y_agr = [0.0, 50.0])
# with a pulse and input rate
par_new = ComponentVector(k_Y = 1e-6, Y0 = 320.0, i_Y = 1e-12,
    ω = 0.01, i_Y_agr = [10.0, 50.0])
# with below ground litter input but negligible decomposition:
# despite lower inputs at the bottom, the stock increases because of downward adv. flux
par_new = ComponentVector(k_Y = 1e-6, Y0 = 220.0, i_Y = 10.0,
    ω = 0.01, i_Y_agr = [0.0, 50.0])
# with decomposition, the pulse barely reaches lower soil
par_new = ComponentVector(k_Y = 1 / 10, Y0 = 80.0, i_Y = 10.0,
    ω = 0.01, i_Y_agr = [0.0, 50.0])
#parl = ComponentVector(prob.p, first(getaxes(par_new)))
#prob2 = remake(prob, p=getdata(par_new), tspan=(0,500))
# reset the initial state
#u0 = Dz_exp.((z_m:dz:-dz), z_m, 2.0) * par_new.Y0
#u0 = Dz_lin.(-(z_m:dz:-dz), -z_m) * par_new.Y0 # straight initial profile
#u0 = Dz_lin.(-(z_m:dz:-dz), -z_m) * par_new.Y0 # straight initial profile
#u0 = Dz_lin.((dz:dz:z_m), z_m) * par_new.Y0 # straight initial profile
#u0 = Dz_lin.((0:dz:z_m)[2:(end-1)], z_m) * par_new.Y0 # straight initial profile
zs = get_1d_grid(prob)
state_pos = get_1d_state_pos(prob)
u0 = ComponentVector(Y = Dz_exp.(zs[state_pos], z_m, 3) * par_new.Y0)
#ax_par = MTKHelpers.axis_of_nums(parameters(get_system(prob)))  

paropt = ComponentVector(state = u0, par = par_new)
pset = ODEProblemParSetter(get_system(prob), paropt)
prob2 = remake(prob, paropt, pset)

#solp = sol = solve(prob2, TRBDF2()); #slower than Rodas5P but more points to plot
#solp = sol = solve(prob2, Rodas5P(), saveat=2); 
#solp = sol = solve(prob2, Rodas5P(), saveat=1, tspan=(70,150)); 
solp = sol = solve(prob2, Rodas5P(); tstops = [80.0]);
solp[t]
solp[Y(t, z)][1, :]
solp[Y(t, z)][end, :]
solp[Yi(t, z)][end, end]

plot(solp[Y(t, z)][end, :], solp[z], xlim = (0, Inf))
plot(solp[Y(t, z)][end, :], solp[z])

# animation of the profile
#(i_t, t_i) = last(enumerate(solp[t]))
Y_ex = extrema(solp[Y(t, z)])
anim = @animate for (i_t, t_i) in enumerate(solp[t])
    plot(solp[Y(t, z)][i_t, :], solp[z], xlim = Y_ex, xlab = "Y(z)", ylab = "z (m)",
        title = "t = $(round(t_i;sigdigits=2))",
        yflip = true)
end;
gif(anim, fps = 3)

# show the input pulse of 50g entering and vanishing again
plot(solp[t], solp[Yi(t, z)][:, end]) # integrated stocks over time

# integration across advection input to soil layer should match agr input
plot(; xlim = (0, 130))
#plot!(solp[t], solp[adv_Yo(t, z)][:,end], label="adv_Y" ) # above ground influx?
#plot!(solp[t], solp[adv_Yo(t, z)][:,end-2] ) 
plot!(solp[t], fagr.(solp[t], par_new.i_Y_agr, par_new.i_Y_agr_pulse), label = "agr input")
#plot!(solp[t], solp[Yi(t, z)][:,end], label="Yi" ) # Y integrated across profile
#plot!(solp[t], -solp[Y_z(t, z)][:,end], label = "Y_z(0)" ) 
plot!(solp[t], solp[adv_Yi(t, z)][:, end], label = "adv_Yi") # Y integrated across profile
#plot!(solp[t], solp[Y(t, z)][:,end]*dz, label="Y_0" ) # Y integrated across profile

i_t = Int(76 / 2)
i_t0 = 1
v = solp[Y(t, z)]
#v = solp[adv_Yo(t,z)]
#v = solp[Y_z(t,z)]
plot(; xlim = extrema(v[(i_t0 + 1):6, :]))
for i in 6:-1:1
    local i_t = i_t0 + i
    plot!(v[i_t, :], solp[z], label = sol[t][i_t])
end
plot!(v[i_t0, :], solp[z], label = sol[t][i_t0])
