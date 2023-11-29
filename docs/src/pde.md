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
get_discrete_space
```

