Pkg.activate(;temp=true)
#Pkg.add(["OrdinaryDiffEq","ModelingToolkit","SymbolicIndexingInterface","ComponentArrays","ForwardDiff"])
Pkg.add(["OrdinaryDiffEq","ComponentArrays","ForwardDiff","SciMLStructures"])
Pkg.add(["Optimization","OptimizationOptimJL"])
Pkg.develop("SymbolicIndexingInterface")
Pkg.develop("ModelingToolkit")

using ModelingToolkit
using ModelingToolkit: t_nounits as t, D_nounits as D

using OrdinaryDiffEq

using SymbolicIndexingInterface: parameter_values, state_values
using SciMLStructures: Tunable, replace, replace!, canonicalize, isscimlstructure

using Optimization
using OptimizationOptimJL


@parameters α β γ δ
@variables x(t) y(t)
eqs = [D(x) ~ (α - β * y) * x
       D(y) ~ (δ * x - γ) * y]
@mtkbuild odesys = System(eqs, t)


odeprob = ODEProblem(
    odesys, [x => 1.0, y => 1.0], (0.0, 10.0), [α => 1.5, β => 1.0, γ => 3.0, δ => 1.0])

@code_warntype replace(Tunable(), parameter_values(odeprob), collect(1.0:4.0))

timesteps = 0.0:0.1:10.0
sol = solve(odeprob, Tsit5(); saveat = timesteps)
data = Array(sol)
# add some random noise
data = data + 0.01 * randn(size(data))

# obtain template and indices of parameters
function make_template_and_pos(odeprob, odesys, par_opt)
       x_template = collect(canonicalize(Tunable(), odeprob.p)[1])
       prob_ind = remake(odeprob, p = replace(Tunable(), 
              parameter_values(odeprob), eachindex(x_template))) 
       pos_opt = prob_ind.ps[par_opt]
       pos_par = zeros(Int,length(x_template))
       for (i,pos) in pairs(pos_opt)
              pos_par[pos] = i
       end
       x_template, pos_par  # zero: indicates non optimized, >0 indicates position in par_opt
end
#par_opt = [α, β, γ, δ]
par_opt = [β, γ, δ]
x_template, pos_par = make_template_and_pos(odeprob, odesys, par_opt)
#par_opt[pos_par]


function loss(popt, p)
       T = eltype(popt)
       odeprob = p[1] # ODEProblem stored as parameters to avoid using global variables
       ps = parameter_values(odeprob) # obtain the parameter object from the problem
       x_template = p[4]
       # construct a new parameters vector from template and optimized parameters
       pos_par = p[5]
       xnew = [iszero(pos_par[i]) ? T(x_template[i]) : popt[pos_par[i]] for i in eachindex(x_template)]::Vector{T}
       # the following is inferred as Any? -> causing newprob, sol, data, sum to Any
       ps = replace(Tunable(), ps, xnew) # create a copy with the values passed to the loss function
       # we also have to convert the `u0` vector
       u0 = T.(state_values(odeprob))
       # remake the problem, passing in our new parameter object
       newprob = remake(odeprob; u0 = u0, p = ps)
       timesteps = p[2]
       sol = solve(newprob, AutoTsit5(Rosenbrock23()); saveat = timesteps)
       truth = p[3]
       data = Array(sol)
       return sum((truth .- data) .^ 2) / length(truth)
   end
   

loss(odeprob.ps[par_opt].*0.8, (odeprob, timesteps, data, x_template, pos_par))      
#using Cthulhu
@descend_code_warntype loss(odeprob.ps[par_opt].*0.8, (odeprob, timesteps, data, x_template, pos_par))   
using ForwardDiff
ForwardDiff.gradient(x -> loss(x, (odeprob, timesteps, data, x_template, pos_par)), odeprob.ps[par_opt].*0.8)


# manually create an OptimizationFunction to ensure usage of `ForwardDiff`, which will
# require changing the types of parameters from `Float64` to `ForwardDiff.Dual`
optfn = OptimizationFunction(loss, Optimization.AutoForwardDiff())
# parameter object is a tuple, to store differently typed objects together
optprob = OptimizationProblem(
       optfn, odeprob.ps[par_opt].*0.8, (odeprob, timesteps, data, x_template, pos_par), lb = 0.1zeros(length(par_opt)), ub = 3ones(length(par_opt)))
sol = solve(optprob, BFGS())
sol.u
odeprob.ps[par_opt]


replace!(Tunable(), parameter_values(odeprob), sol.u)
odeprob.ps[[α, β, γ, δ]]



canonicalize(Tunable(), odeprob.p)[1]

