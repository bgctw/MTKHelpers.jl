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
