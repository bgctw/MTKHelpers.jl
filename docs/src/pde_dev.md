```@meta
CurrentModule = MTKHelpers
```

# PDE support developer notes

## Depth distributions
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
\\
F(x) &= \frac{1}{e^{-bx_m} - 1} \left( e^{-bx} -1 \right)
\end{aligned}
```

For functions and figures see [Exponential depth distribution of inputs](@ref).

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

For positive direction z we have

```math
\begin{align}
-\omega \frac{\partial y_0}{\partial z} 
&\approx - \frac{\omega y_0 - \omega y_{-1}}{\Delta z} 
\\
&= - \frac{\omega y_0 - f(t)}{\Delta z} 
\\
&= \frac{f(t) - \omega y_0}{\Delta z} 
\\
\frac{\partial y_0}{\partial z}  &\approx \frac{y_0 - f(t)/\omega}{\Delta z} 
\end{align}
```


Note that ``y(z)`` is a density per unit ``z`` (here ``g``):  ``[y] = g/m``. Hence 
``\omega y_{n+1}`` has units ``\frac{m}{s} \frac{g}{m} = \frac{g}{s}``, i.e. th
units of a flux in time.
