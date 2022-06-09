using ModelingToolkit, DifferentialEquations
@variables t 
D = Differential(t) 
sts = @variables x(t) RHS(t)  # RHS is observed
ps = @parameters τ=0.3      # parameters
@named m = ODESystem([ RHS  ~ (1 - x)/τ, D(x) ~ RHS ], t, sts, ps)
ms = structural_simplify(m)
prob = ODEProblem(ms, [x => 1.1], (0.0,1.0), [])
sol = solve(prob, Tsit5())

   
using BenchmarkTools
const constant = 10
# actually constant can be reassigned by the same type, only issues warnings
add(x) = constant + x
@btime add(3)

var = 10
add2(x) = var::Int + x
@btime add2(3)

add3(x, y=var) = var +x
@btime add3(3, $var)

const semi = Ref(10)
add4(x) = x + semi[]
@btime add4(3)


using ModelingToolkit, DifferentialEquations, LabelledArrays
using MTKHelpers
function samplesystem(;name,τ = 3.0, p1=1.1, p2=1.2) 
    @variables t 
    D = Differential(t) 
    sts = @variables x(t) RHS(t)             # RHS is observed
    ps = @parameters τ=τ p1=p1 p2=p2       # parameters
    ODESystem([ RHS  ~ p1/p2 * (1 - x)/τ, D(x) ~ RHS ], t, sts, ps; name)
end                     
@named m = samplesystem()
@named sys = embed_system(m)
prob = ODEProblem(sys, [m.x => 0.0], (0.0,10.0))

# setup position matching#
popt = SLVector(m₊x=0.1, m₊p1=2.1)
ps = ProblemParSetter(sys, keys(popt))

# extract optimized 
get_paropt_labeled(ps, prob)

# update states and parameters
prob2 = update_statepar(ps, popt, prob)
get_paropt_labeled(ps, prob2) == popt

# Documentation example
using ModelingToolkit, DifferentialEquations, ComponentArrays
using MTKHelpers
function samplesystem(;name,τ = 3.0, p1=1.1, p2=1.2) 
    @variables t 
    D = Differential(t) 
    sts = @variables x(t) RHS(t)             # RHS is observed
    ps = @parameters τ=τ p1=p1 p2=p2       # parameters
    ODESystem([ RHS  ~ p1/p2 * (1 - x)/τ, D(x) ~ RHS ], t, sts, ps; name)
end                     
@named m = samplesystem()
@named sys = embed_system(m)
prob = ODEProblem(sys, [m.x => 0.0], (0.0,10.0))

# setup position matching
popt = ComponentVector(m₊x=0.1, m₊p1=2.1)
#ps = ProblemParSetter(sys, first(getaxes(popt))) # better type-stability
ps = ProblemParSetter(sys, popt)

# extract optimized 
get_paropt_labeled(ps, prob)

# update states and parameters
prob2 = update_statepar(ps, popt, prob)
prob2.p # p is still a plain vector
label_par(ps, prob2.p).m₊p1 == popt.m₊p1 # attach labels and access properties
label_state(ps, prob2.u0).m₊x == popt.m₊x # attach labels and access properties
get_paropt_labeled(ps, prob2) == popt
