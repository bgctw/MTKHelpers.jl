using ModelingToolkit
#using DifferentialEquations
using OrdinaryDiffEq
#using Plots
using SymbolicIndexingInterface
using SciMLStructures: replace!, Tunable, canonicalize

function lotka_volterra(;name)
    @parameters t
    pars = @parameters α β γ δ
    D = Differential(t)
    vars = @variables x(t) y(t)
    eqs = [
        D(x) ~ (α - β * y) * x, # prey
        D(y) ~ (δ * x - γ) * y # predator
    ]
    tspan = (0.0, 10.0)
    sys = System(eqs, t, vars, pars; tspan = tspan, name = name)
    return sys
end

@mtkbuild lvsys = lotka_volterra();
prob = ODEProblem(lvsys, [lvsys.x => 1.0, lvsys.y => 1.0], lvsys.tspan, [lvsys.α => 1.5, lvsys.β => 1.0, lvsys.γ => 3.0, lvsys.δ => 1.0]);

# Plot simulation.
# plot(solve(prob, Tsit5()))
sol = solve(prob, Tsit5(); saveat=0.1);
data = Array(sol) + 0.01 * randn(size(Array(sol)));

function loss(x, p)
    prob = p[1]
    prob = remake(prob, p = [:α => x[1], :β => x[2], :γ => x[3], :δ => x[4]])
    sol = solve(prob, Tsit5(), saveat=0.1);
    return sum((sol .- data).^2)
end

ForwardDiff.gradient(x -> loss(x,[prob]), rand(4))

using Optimization
using OptimizationOptimJL

f = OptimizationFunction(loss, Optimization.AutoForwardDiff())
oprob = OptimizationProblem(f, rand(4), [prob]);
sol = solve(oprob, BFGS())