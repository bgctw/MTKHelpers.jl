using ModelingToolkit, DifferentialEquations
@variables t 
D = Differential(t) 
sts = @variables x(t) RHS(t)  # RHS is observed
ps = @parameters τ=0.3      # parameters
@named m = ODESystem([ RHS  ~ (1 - x)/τ, D(x) ~ RHS ], t, sts, ps)
ms = structural_simplify(m)
prob = ODEProblem(ms, [x => 1.1], (0.0,1.0), [])
sol = solve(prob)

sol(0.3, idxs=[RHS])                                                                                                                                                                     
