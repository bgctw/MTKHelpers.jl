function samplesystem(;name) 
    @variables t 
    D = Differential(t) 
    sts = @variables x(t) RHS(t)  # RHS is observed
    ps = @parameters τ       # parameters
    ODESystem([ RHS  ~ (1 - x)/τ, D(x) ~ RHS ], t, sts, ps; name)
end                     

@testset "embed with name different than m" begin
    @named m2 = samplesystem()
    @named sys = embed_system(m2)
    prob = ODEProblem(sys, [m2.x => 0.0], (0.0,10.0), [m2.τ => 3.0])
    sol = solve(prob)
    @test first(sol[m2.x]) == 0.0
    #plot(sol, vars=[m2.x,m2.RHS])    
end
