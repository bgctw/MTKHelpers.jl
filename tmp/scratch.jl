using MTKHelpers


using ModelingToolkit, DifferentialEquations, StaticArrays
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

popt = SLVector(m₊x=0.1, m₊p1=2.1)
ps = ProblemParSetter(sys, keys(popt))

get_paropt(ps, prob; label=true)

prob2 = update_statepar(ps, popt, prob)
get_paropt(ps, prob2; label=true) == popt

prob2.u0[1] == popt.m₊x
prob2.p[2] == popt.m₊p1

