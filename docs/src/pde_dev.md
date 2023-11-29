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

