# develop new ODEParSetter that is based on populating canonical buiffers and replacing portions

import Pkg; 
Pkg.add(["OrdinaryDiffEq","Optimization","OptimizationPolyalgorithms","SciMLSensitivity",
  "Zygote","ForwardDiff","Plots","ModelingToolkit","SymbolicIndexingInterface",
  "SciMLStructures","SciMLBase","Plots"])
import OrdinaryDiffEq as ODE
import Optimization as OPT
import OptimizationPolyalgorithms as OPA
import SciMLSensitivity as SMS
import Zygote
import ForwardDiff
using SciMLBase
#@usingany Plots: Plots
#using Plots: Plots

using ModelingToolkit, OrdinaryDiffEq
using ModelingToolkit: t_nounits as t, D_nounits as D, ModelingToolkit as MTK
using SymbolicIndexingInterface: SymbolicIndexingInterface as SII
using SciMLStructures: SciMLStructures as SS

@variables x(t) y(t) z(t)
@parameters p_lv[1:2]=[1.5,1.0] γ=3.0 δ=1.0

eqs = [D(x) ~ p_lv[1] * x - p_lv[2] * x * y
       D(y) ~ -γ * y + δ * x * y
       z ~ x + y]

@named sys = System(eqs, t)
simpsys = mtkcompile(sys)
tsteps = 0.0:0.1:10.0
tspan = extrema(tsteps)
prob = ODEProblem(simpsys, [x => 1.1, y => 1.2], tspan)
#prob = ODEProblem(simpsys, [], tspan; guesses=[x => 1.1, y => 1.2])

sol = sol_true = ODE.solve(prob, ODE.Tsit5(), saveat = tsteps)
# Plot the solution
#Plots.plot(sol)

#paropt = [x, β]  # currently initial state x cannot be handled by SII.setp
paropt = [p_lv,γ]
stateopt = [x]
p_true = vcat(prob.ps[paropt]..., prob[stateopt]...) # also optimize initial state x
p0 = p = p_true .+ randn(length(p_true)) .* 0.2
probo = remake(prob)


ps_ind = MTK.parameter_index.(Ref(simpsys), paropt)
x_ind = MTK.variable_index.(Ref(simpsys), stateopt) # 2

#--------------- SciMLStructures
# TODO: describe how post-modify tunable metainformation in the system
#    to not rely on the user constructing the buf and bufx objects
#    and buffers only pertain to really tunable parameters
pv = SII.parameter_values(prob)
buf, _, _ = SS.canonicalize(SS.Tunable(), pv)
bufx, _, _ = SS.canonicalize(SS.Initials(), pv)


par_all = parameters(simpsys)
ps_ind = Dict(par_all .=> MTK.parameter_index.(Ref(simpsys), par_all))
# can assume that indices are ordered
#buf_starts = [minimum(ps_ind[k].idx) for k in par_all]

# todo get symbolic sys.pvar by symbol psym
par:opt 

ends = cumsum(length.(paropt))
starts = vcat(1, ends[1:(end-1)] .+ 1)
ppos_opt = Dict(paropt[i] => starts[i]:ends[i] for i in axes(paropt,1))
# construct Tuples with entries for each par_all
buf_ind = Tuple((MTK.parameter_index.(Ref(simpsys), par_all)))
p_ind = Tuple((get(ppos_opt, k, missing) for k in par_all))
is_tuned = par_all .∈ Ref(keys(ppos))

pv = SII.parameter_values(prob)
buf, _, _ = SS.canonicalize(SS.Tunable(), pv)
bufx, _, _ = SS.canonicalize(SS.Initials(), pv)
ET = promote_type(eltype(p), eltype(buf))
bufN_v = (is_tuned[i] ? p[p_ind[i]] : buf[buf_ind[i].idx] for (i,k) in enumerate(par_all))
bufn = vcat(bufN_v...)::Vector{ET}
pvn = SS.replace(SS.Tunable(), pv, bufn)
#
# care for initials
probo = remake(probo; p = pvn)



k = p_lv


getindex(ppos,δ, nothing)



function compute_sol(p)
    # local pv = SII.parameter_values(probo)
    # local buf, _, _ = SS.canonicalize(SS.Tunable(), pv)
    # local bufx, _, _ = SS.canonicalize(SS.Initials(), pv)
    local bufn = vcat(p[1:2], buf[3:4]) 
    local pvn = SS.replace(SS.Tunable(), pv, bufn)
    #local pvn2 = pvn
    local xvn = vcat(xv[1], p[3])
    local bufxn = vcat(bufx[1], p[3], bufx[3:end]) # index 2?
    local pvn2 = SS.replace(SS.Initials(), pvn, bufxn)
    #local probon = remake(probo; u0 = xvn, p = pvn) # u0 set to pvn.initials instead
    local probon = remake(probo; p = pvn2) 
    local probon2 = probon
    #local probon2 = remake(probon; u0 = xvn) # mutating?
    #@show probon2.u0, probon2.p
    local sol = ODE.solve(probon2, ODE.Tsit5(), saveat = tsteps)
    return sol
end

function loss(p)
    local sol = compute_sol(p)
    local loss = SciMLBase.successful_retcode(sol) ? sum(abs2, sol .- sol_true) : 1e30
    return loss
end

#include("tmp/test.jl")
loss(p0)
loss(p)
loss(p_true)

callback = function (state, l)
    display(l)
    p = state.u
    sol = compute_sol(p)
    plt = Plots.plot(sol, ylim = (0, 7))
    display(plt)
    # Tell Optimization.solve to not halt the optimization. If return true, then
    # optimization stops.
    return false
end

adtype = OPT.AutoForwardDiff()
adtype = OPT.AutoZygote()
optf = OPT.OptimizationFunction((x, p) -> loss(x), adtype)
#optf = OPT.OptimizationFunction((x, p) -> loss(x))
optprob = OPT.OptimizationProblem(optf, p0)

opt_alg = OPA.PolyOpt()
#opt_alg = OPA.LBFGS()
#include("tmp/test.jl")
result_ode = OPT.solve(optprob, opt_alg,
    callback = callback,
    maxiters = 10,
    )

result_ode.stats
p = result_ode.u
hcat(p0, p, p_true)
loss(result_ode.u)


