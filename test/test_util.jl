@testset "symbols_state ODE" begin
    @named m = samplesystem()
    @test symbols_state(m) == [:x, :RHS]   
    @test symbols_par(m) == [:τ, :p1, :p2]  
end

@testset "symbol" begin
    @named m = samplesystem()
    @test symbol(m.x) == :m₊x # calls Num which calls Par
    @test symbol(m.τ) == :m₊τ 
end;

@testset "strip_namespace" begin
    @test strip_namespace(:x) == :x
    @test strip_namespace(:m₊x) == :x
    @test strip_namespace(:s₊m₊x) == :x
    @test strip_namespace("s.m.x") == "x"
end;


@testset "embed with name different than m" begin
    @named m2 = samplesystem()
    @named sys = embed_system(m2)
    prob = ODEProblem(sys, [m2.x => 0.0], (0.0,10.0), [m2.τ => 3.0])
    sol = solve(prob)
    @test first(sol[m2.x]) == 0.0
    #plot(sol, vars=[m2.x,m2.RHS])    
end

