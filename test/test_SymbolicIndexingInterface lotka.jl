import Pkg; 
Pkg.activate(;temp=true)
Pkg.add(["OrdinaryDiffEq","Optimization","OptimizationPolyalgorithms","SciMLSensitivity",
  "Zygote","ForwardDiff","ModelingToolkit","SymbolicIndexingInterface",
  "SciMLStructures","SciMLBase","ComponentArrays","Plots"])
import OrdinaryDiffEq as ODE
import Optimization as OPT
import OptimizationPolyalgorithms as OPA
import SciMLSensitivity as SMS
import Zygote
import ForwardDiff
using SciMLBase
using ComponentArrays
using ComponentArrays: ComponentArrays as CA
#@usingany Plots: Plots
#using Plots: Plots

using ModelingToolkit, OrdinaryDiffEq
using ModelingToolkit: t_nounits as t, D_nounits as D, ModelingToolkit as MTK
using SymbolicIndexingInterface: SymbolicIndexingInterface as SII
using SciMLStructures: SciMLStructures as SS

@variables x(t) y(t) z(t)
@parameters px[1:2]=[1.5, 1.0] γ=3.0 δ=1.0

eqs = [D(x) ~ px[1] * x - px[2] * x * y
       D(y) ~ -γ * y + δ * x * y
       z ~ x + y]

@named sys = System(eqs, t)
simpsys = mtkcompile(sys)
tsteps = 0.0:0.1:10.0
tspan = extrema(tsteps)
prob = ODEProblem(simpsys, [x => 1.1, y => 1.2], tspan)
sol = sol_true = ODE.solve(prob, ODE.Tsit5(), saveat = tsteps)
#Plots.plot(sol)

paropt = [px, γ]
paropt_sym = Symbol.(paropt)
n_paropt = sum(length.(paropt))
stateopt = [x]
stateopt_sym = Symbol.(stateopt)
#p_true = vcat(prob.ps[paropt]..., prob[stateopt]...) 
p_true = vcat(
    ComponentVector(;zip(paropt_sym,prob.ps[paropt])...),
    ComponentVector(;zip(stateopt_sym,prob[stateopt])...))
p0 = p = p_true .+ randn(length(p_true)) .* 0.2
probo = remake(prob) # copy
nest_structure(p, syms) = [p[k] for k in syms]
parstateopt = vcat(paropt, stateopt)
parstateopt_sym = Symbol.(parstateopt)
p0n = nest_structure(p0, parstateopt_sym)

#--------------- recommended method for optimization at MTK.FAQ
# obstacle: do not find documentation for initial conditions
# obstacle: problems with both Zygote and ForwardDiff
setter! = SII.setp(simpsys, paropt)
setter!(probo, nest_structure(p0, paropt_sym)) 
function loss(p)
    local p_struc = nest_structure(p, paropt_sym) #  omit u0
    setter!(probo, p_struc)
    local sol = ODE.solve(probo, ODE.Tsit5(), saveat = tsteps)
    local loss = sum(abs2, sol .- sol_true)
    return loss
end
loss(p0)
#Zygote.gradient(loss, p0)
#ForwardDiff.gradient(loss, p0)


#------------- alternative using indices, 
# documentation for initial values
# obstacle: problems with both Zygote and ForwardDiff 
probo = remake(prob)
ps_ind = MTK.parameter_index.(Ref(simpsys), paropt)
setindex!.(Ref(probo.ps), nest_structure(p0, paropt_sym), ps_ind)
probo.ps[px] == p0[Symbol("px[1:2]")]
x_ind = MTK.variable_index.(Ref(simpsys), stateopt) # 2
#setindex!(probo.ps, [p0[3]], x_ind)

function loss(p)
    setindex!.(Ref(probo.ps), nest_structure(p, paropt_sym), ps_ind)
    # for simplicity, start omitting u0
    local xv = probo.u0
    # local xvn = [begin
    #     ii = findfirst(==(i), x_ind)
    #     isnothing(ii) ? xv[i] : p[length(ps_ind) + ii]
    # end for i in  axes(xv,1)]
    #local xvn = vcat(xv[1], p[3])
    #probox = remake(probo, u0 = xvn) 
    local sol = ODE.solve(probo, ODE.Tsit5(), saveat = tsteps)
    local loss = sum(abs2, sol .- sol_true)
    return loss
end
loss(p0)
Zygote.gradient(loss, p0) # nothing?
#ForwardDiff.gradient(loss, p0)


#--------------- alternative: SciMLStructures
# works with ForwardDiff
# obstacle: how to infere positions of optimized parameters in canicalized buffer?
# obstacle: error or zero Zygote.gradient for initial conditions ?
probo = remake(prob)
function loss(p)
    local sol = compute_sol(p, probo)
    local loss = SciMLBase.successful_retcode(sol) ? sum(abs2, sol .- sol_true) : 1e30
    return loss
end
xv = probo.u0
function compute_sol(p, probo) # variant without Zygote compile error, but not general
    local pv = SII.parameter_values(probo)
    local buf, _, _ = SS.canonicalize(SS.Tunable(), pv)
    local bufx, _, _ = SS.canonicalize(SS.Initials(), pv)
    local bufn = vcat(p[1:2], buf[3], p[3])       # TODO describe general
    local pvn = SS.replace(SS.Tunable(), pv, bufn)
    local bufxn = vcat(bufx[1], p[4], bufx[3:end]) # TODO describe general
    local pvn2 = SS.replace(SS.Initials(), pvn, bufxn)
    #local probon = remake(probo; u0 = xvn, p = pvn) # u0 set to pvn.initials instead
    local probon = remake(probo; p = pvn2) 
    local probon2 = probon
    #local xvn = vcat(xv[1], p[3])
    #local probon2 = remake(probon; u0 = xvn) # mutating?
    #@show probon2.u0, probon2.p
    local sol = ODE.solve(probon2, ODE.Tsit5(), saveat = tsteps)
    return sol
end
psetter! = SII.setp(probo, paropt)
ssetter! = SII.setu(probo, stateopt)
function compute_sol(p, probo)
    local pv = SII.parameter_values(probo)
    local buf, _, _ = SS.canonicalize(SS.Tunable(), pv)
    local bufx, _, _ = SS.canonicalize(SS.Initials(), pv)
    # for ForwardDiff need to convert the eltype of the entire portions
    ET = eltype(p)
    local pvt = pv
    pvt = eltype(buf) == ET ? pvt : SS.replace(SS.Tunable(), pvt, ET.(buf)) 
    pvt = eltype(bufx) == ET ? pvt : SS.replace(SS.Initials(), pvt, ET.(bufx)) 
    #local psetter! = SII.setp(probo, paropt)
    local p_struc = nest_structure(p, paropt_sym) 
    psetter!(pvt, p_struc)
    #local ssetter! = SII.setu(probo, stateopt)
    local s_struc = nest_structure(p, stateopt_sym) 
    #ssetter!(pvt, s_struc)   # state_values(MTKParameters) not implemented
    ssetter!(probon2, s_struc)
    #local sol = ODE.solve(probon2, ODE.Tsit5(), saveat = tsteps)
    # need another remake to update probon2.p.initial to probon2.u0
    local probon3 = remake(probon2, u0=probon2.u0)
    #@show probon3.u0, probon3.p
    local sol = ODE.solve(probon3, ODE.Tsit5(), saveat = tsteps)
    return sol
end
#include("tmp/test.jl")
loss(p0)
ForwardDiff.gradient(loss, p0)
#Zygote.gradient(loss, p0)
#Zygote.gradient(loss, CA.getdata(p0))
loss(p)
loss(p_true)

callback = function (state, l)
    display(l)
    # p = state.u
    # sol = compute_sol(p)
    # plt = Plots.plot(sol, ylim = (0, 7))
    # display(plt)
    # Tell Optimization.solve to not halt the optimization. If return true, then
    # optimization stops.
    return false
end

adtype = OPT.AutoForwardDiff()
#adtype = OPT.AutoZygote()
optf = OPT.OptimizationFunction((x, p) -> loss(x), adtype)
optprob = OPT.OptimizationProblem(optf, p0)

opt_alg = OPA.PolyOpt()
#opt_alg = OPA.LBFGS()
result_ode = OPT.solve(optprob, opt_alg,
    callback = callback,
    maxiters = 20,
    )

result_ode.stats
p = result_ode.u
hcat(p0, p, p_true)
loss(result_ode.u)
