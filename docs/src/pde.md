```@meta
CurrentModule = MTKHelpers
```

# PDE suport

## Exponential depth distribution

The inputs to the system have to be described by a differential across depth.
The following derivation provides an exponential function ``f(z;b) = a e^{bz}``,
such that ``\int_{z_m}^0 f(z) dz = 1``.
Such a function can be used to distribute total inputs across depth.

Derivation: for positive ``x = -z``:

```math
\begin{aligned}
f(x) &= a e^{-bx}
\\
F(x) &= \frac{a}{-b} e^{-bx} 
\\
\left[ F(x) \right]_0^{x_m} &= \frac{a}{-b} \left( e^{-b x_m} -1 \right) = 1
\\
a &= \frac{b}{1 - e^{-bx_m}}
\\
f(x) &= \frac{b}{1 - e^{-bx_m}} e^{-bx}
\end{aligned}
```

Hence 
```math
\begin{aligned}
f(z) &= \frac{b}{1 - e^{b z_m}} e^{b z}
\\
F(z) &= \frac{1}{1 - e^{b z_m}} e^{b z}
\end{aligned}
```

The following figure displays the input at depth z
```@example Dz_exp
using Plots
z_m = -2
dz = 0.1
z = z_m:dz:0
Dz_exp(z,z_m,b) = b/(1-exp(b*z_m)) * exp(b*z)
# Integral from z_m to z
Iz_exp(z,z_m,b) = 1/(1-exp(b*z_m)) * (exp(b*z)-exp(b*z_m))
Iz_exp(0.0, z_m, 0.5) == 1
Dz_exp.(z, z_m, 0.5)
Iz_exp.(z, z_m, 0.5) 
# Check by approximating the integral by sum
Iz_emp(z, z_m, args...; dz=1e-3) = sum(z_i -> Dz_exp(z_i,z_m,args...)*dz, z_m:dz:z)
isapprox(Iz_emp(0, z_m, 2), 1; atol=1e-2)
plot(Dz_exp.(z, z_m, 0.5), z, ylab="depth z (m)", xlab="di/dz", label="b = 0.5")
plot!(Dz_exp.(z, z_m, 1), z, label="b=1.0")
plot!(Dz_exp.(z, z_m, 2), z, label="b=2.0")
```

And the following figure the integral over the inputs from ``z_m`` to ``z``.
```@example Dz_exp
plot(Iz_exp.(z, z_m, 0.5), z, ylab="depth z (m)", xlab="integral input from z_m to z", label="b = 0.5")
plot!(Iz_exp.(z, z_m, 2), z, label="b=2.0")
```

## Input at upper boundary

The advetion term is the same when changing direction e.g. ``x \in [0 \ldots x_m]`` 
with velocity ``\omega_x > 0``
versus ``z \in [z_m \ldots 0]`` with velocity ``\omega = -\omega_x < 0``.
Note that also ``\Delta y = y_{i+1} - y_{i} = -\Delta x`` and ``\partial y = -\partial x``.

```math
\begin{align}
\frac{\partial y}{\partial t} + \omega_x \frac{\partial y}{\partial x} = 0
\\
\operatorname{adv} 
=  -\omega_x \frac{\partial y}{\partial x} 
=  -\omega \frac{\partial y}{\partial z} 
= \frac{\partial y}{\partial t}
\end{align}
```

The advection term is discretized with the upwind scheme using the state of the
layer that comes before the current state. For ``\omega < 0`` this is state at ``i+1``.
```math
\operatorname{adv} 
\approx -\omega \frac{y_{i+1} - y_i}{\Delta z} 
```

The idea for a boundary condition that describes
a fixed flux into the system at ``z_n = 0`` is to replace
the flux from the upwind state by the influx ``-f(t)``. 
The minus sign indicates that the direction goes opposite to positive ``z``.

```math
\begin{align}
-\omega \frac{\partial y_n}{\partial z} 
&\approx - \frac{\omega y_{n+1} - \omega y_n}{\Delta z} 
\\
&= - \frac{-f(t) - \omega y_n}{\Delta z} 
\\
&= \frac{f(t) + \omega y_n}{\Delta z} 
\\
\frac{\partial y_n}{\partial z}  &\approx \frac{-f(t)/\omega - y_n}{\Delta z} 
\end{align}
```

Note that ``y(z)`` is a density per unit ``z`` (here ``g``):  ``[y] = g/m``. Hence 
``\omega y_{n+1}`` has units ``\frac{m}{s} \frac{g}{m} = \frac{g}{s}``, i.e. th
units of a flux in time.
