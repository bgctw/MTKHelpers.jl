```@meta
CurrentModule = MTKHelpers
```

# PDE support

## Exponentially increasing grid

```@docs
grid_exp
```

```@example doc
using MTKHelpers, Plots
n = 16
z_m = 0.3
scatter(grid_exp(n, z_m, 2), label="efold=2");
scatter!(grid_exp(n, z_m, 4), label="efold=4") 
```


## Exponential depth distribution of inputs

The inputs to the 1D PDE system have to be described by an unnormalized density.
The following function decreases exponentially from zero to ``x = -z`` with
with an integral of 1. It can be used to distribute the total input across the domain.

```@docs
Dz_exp
Dz_lin
```

For its derivation see [PDE support developer notes](@ref).

The following figure displays the density for several e-folding times.

```@example Dz_exp
using MTKHelpers, Plots
x_m = 2.0
dx = 0.1
x_grid = 0:dx:x_m
Iz_exp(x_m, x_m, 0.5) == 1
# Check by approximating the integral by sum
Iz_emp(x, x_m, args...; dx=1e-3) = sum(x_i -> Dz_exp(x_i,x_m,args...)*dx, 0:dx:x)
isapprox(Iz_emp(x_m, x_m, 2), 1; atol=1e-2)
plot(;ylab="depth x (m)", xlab="∂/∂x", yflip=true);
plot!(Dz_exp.(x_grid, x_m, 0.5), x_grid, label="b = 0.5");
plot!(Dz_exp.(x_grid, x_m, 1), x_grid, label="b=1.0");
plot!(Dz_exp.(x_grid, x_m, 2), x_grid, label="b=2.0")
```

And the following figure shows the corresponding integral.

```@example Dz_exp
plot(;ylab="depth x (m)", xlab="integral input from 0 to x_m", yflip=true);
plot!(Iz_exp.(x_grid, x_m, 0.5), x_grid, label="b = 0.5");
plot!(Iz_exp.(x_grid, x_m, 2), x_grid, label="b=2.0")
```

For negative z, these functions can be used with the minus sign

```@example Dz_exp
z_m = -2
dz = 0.1
z_grid = z_m:dz:0
Iz_exp(-z_m, -z_m, 0.5) == 1
plot(;ylab="depth z (m)", xlab="di/dz");
plot!(Dz_exp.(-z_grid, -z_m, 0.5), z_grid, label="b = 0.5");
plot!(Dz_exp.(-z_grid, -z_m, 1), z_grid, label="b=1.0");
plot!(Dz_exp.(-z_grid, -z_m, 2), z_grid, label="b=2.0")
```
## Extracting grid information from the system or problem

The following require MethodOfLines to be loaded to activate the extension.
```@docs
get_1d_grid
get_1d_state_pos
get_discrete_space
```

## Example pde

Define the system, here with using a non-homogeneous grid ([`grid_exp`](@ref)) 
and using symbolic array parameters.

```@example pde
using OrdinaryDiffEq, ModelingToolkit, DomainSets
using MethodOfLines
using MTKHelpers
using LoggingExtras # hide

@parameters t z 
z_m = +0.3 # maximum depth in m, change to positive depth
n_z = 16
#z_grid = collect(range(0,z_m, length=n_z))
z_grid = grid_exp(n_z, z_m, 3.0) 
dz = z_m/(n_z-1) # control of same grid represented by a single number
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
    adv_Yi(t,z) ~ Iz(adv_Yo(t, z)), # integral over advection inputs
]
```
The below-ground inputs decrease exponentially with depth using [`Dz_exp`](@ref).

The above-ground input is specified as a function with a baseline and
a pulse at t=80 and connected to the model by a Dirichlet boundary conditions.

The initial values are specified constant across depth using [`Dz_lin`](@ref).

```@example pde; output=false
i_Yz(t,z,i_Y) = Dz_exp(z, z_m,4.0) * i_Y  
#@register_symbolic i_Yz(t, z, i_Y)

eqs = vcat(eqs0, [
    i_Yo(t, z) ~ i_Yz(t, z, i_Y),         # observable input Y
])

fgauss(x, μ, σ2) = 1/sqrt(2*pi*σ2) * exp(-(x-μ)^2/(2*σ2))
fagr(t, i_Y_agr, i_Y_agr_pulse) = i_Y_agr + i_Y_agr_pulse*fgauss(t, 80, 2^2)
#@register_symbolic fagr(t, i_Y_agr, i_Y_agr_pulse)

bcs = [
    Y(0, z) ~ Dz_lin(z, z_m) * Y0, # initial constant with depth,
    ω * Y(t, 0) ~ fagr(t,i_Y_agr[1], i_Y_agr[2]), # specified flux at upper boundary
    ∂_z(Y(t, z_m)) ~ 0, # negligible change in concentration at lower boundary
    #
    # following are only necessary because observables boundaries 
    # need to be specified with vector grid
    ∂_z(dec_Y(t, 0)) ~ 0.0, 
    ∂_z(dec_Y(t, z_m)) ~ 0.0, 
    ∂_z(i_Yo(t, 0)) ~ 0.0, 
    ∂_z(i_Yo(t, z_m)) ~ 0.0, 
    #∂_z(adv_Yo(t, 0)) ~ 0.0, 
    ω * adv_Yo(t, 0) ~ fagr(t,i_Y_agr[1], i_Y_agr[2]), 
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
domains = [t ∈ Interval(0.0, 500.0), z ∈ Interval(0.0, z_m), ]

# PDE system
state_vars = [Y(t, z), i_Yo(t,z), adv_Yo(t,z), dec_Y(t,z), Y_t(t,z), Yi(t,z), i_Yi(t,z), Y_z(t,z), adv_Yi(t,z) ]
params = [
        k_Y => 2.0,
        Y0 => 200.0,
        i_Y => 50.0,
        ω => 0.01,
        i_Y_agr[1] => 10.0, # base input
        i_Y_agr[2] => 50.0, # pulse at t=80       
]
pdesys, prob = LoggingExtras.withlevel(Logging.Error) do 
    @named pdesys = PDESystem(eqs, bcs, domains, [t, z], state_vars, params)
    prob = discretize(pdesys, discretization) 
    pdesys, prob
end; 
nothing # hide
```

Inside a cost function one wants to modify the parameters and initial states
of the problem without the need to reconstruct the problem from the system.

New parameters can be conveniently specified by a ComponentVector of only
those components that need to be updated.

The new initial states must be computed at the grid positions of states, i.e.
at grid positions that are not at the boundary. Here we again
specify a constant profile using [`Dz_lin`](@ref).

```@example pde
using ComponentArrays
# ω not specified
par_new = ComponentVector(k_Y = 1/10, Y0 = 80.0, i_Y = 10.0,  i_Y_agr = [2.0, 50.0]) 
#
zs = get_1d_grid(prob)
state_pos = get_1d_state_pos(prob)
u0 = ComponentVector(Y = Dz_lin.(zs[state_pos], z_m) * par_new.Y0)

paropt = ComponentVector(state = u0, par = par_new)
pset = ODEProblemParSetter(get_system(prob), paropt)
prob2 = remake(prob, paropt, pset)
nothing # hide
```

We can check the new parameters

```@example pde
label_par(pset, prob2.p)[keys(par_new)] == par_new
label_state(pset, prob2.u0)[keys(u0)] == u0
```

And solve the problem and display the results.
```@example pde
sol = solve(prob2, Rodas5P(); tstops=[80.0]); 
nothing # hide
```

```@example pde
using Plots # hide
Y_ex = extrema(sol[Y(t,z)]) # hide
anim = @animate for (i_t, t_i) in enumerate(sol[t]) # hide
    Plots.plot(sol[Y(t,z)][i_t,:], sol[z], xlim=Y_ex, xlab="Y(z)", ylab="z (m)",  # hide
    title="t = $(round(t_i;sigdigits=2))", # hide
    yflip=true, # hide
    ) # hide
end; # hide
gif(anim, fps=3) # hide
```
