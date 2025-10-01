#TestEnv.activate()
using Test
using OrdinaryDiffEq, ModelingToolkit
using ComponentArrays: ComponentArrays as CA
#using StaticArrays: StaticArrays as SA

using SymbolicIndexingInterface: SymbolicIndexingInterface as SII

using ForwardDiff: ForwardDiff
import Zygote
using Optimization
import SciMLSensitivity as SMS

pkgdir = dirname(dirname(pathof(MTKHelpers)))
include(joinpath(pkgdir,"test","samplesystem.jl"))
include(joinpath(pkgdir,"test","testset_utils.jl"))

@named m = samplesystem()

(u1, p1, popt1s, prob_sys1) = get_sys_ex_scalar();
popt1 = flatten1(popt1s)
ps1ps = ps1 = ODEProblemParSetter(get_system(prob_sys1), popt1)
get_paropt_labeled(ps1, prob_sys1)

    @named me = embed_system(m)
    ts = 0.0:0.1:1.0
    prob = ODEProblem(me, vcat([m.x => 1.1]), extrema(ts), saveat = ts, jac=true)
    sol_true = solve(prob, Tsit5())
    #@usingany Plots
    #plot(sol_true[m.x])

paropt = [m.x, m.p2]
paropt = [m.p1, m.p2]
setter! = SII.setp(me, paropt)
p_true = [prob[m.x], prob.ps[m.p2]]
p_true = prob.ps[paropt]
p0 = p_true .+ randn(length(p_true)) .* 0.05
probo = remake(prob)

function loss(p)
    setter!(probo, p)
    sol = solve(probo, Tsit5())
    sum(abs2, sol .- sol_true)
end    
#plot!(sol[m.x])
setter!(probo, p0)
probo.ps

callback = function (state, l)
    display(l)
    pred = ODE.solve(prob, ODE.Tsit5(), p = state.u, saveat = tsteps)
    plt = Plots.plot(pred, ylim = (0, 6))
    display(plt)
    # Tell Optimization.solve to not halt the optimization. If return true, then
    # optimization stops.
    return false
end
optf = OptimizationFunction((x,p) -> loss(x), AutoZygote())
optf = OptimizationFunction((x,p) -> loss(x), AutoForwardDiff())
opt_prob = OptimizationProblem(optf, p0)
opt_sol = solve(opt_prob, Optimization.LBFGS(), maxiters = 100)