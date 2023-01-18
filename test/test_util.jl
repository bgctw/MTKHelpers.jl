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
end;

@testset "override_system" begin
    @named mc = samplesystem_const(-0.1)
    @named sys = embed_system(mc)
    prob = ODEProblem(sys, [mc.x => 1.0], (0.0,10.0))
    sol = solve(prob)
    #@test sol(8)[1] ≈ 1-0.1*8
    @test isapprox(sol(8)[1], exp(-0.1*8), atol = 1e-5)
end;

@testset "override_system: error on unknown equation" begin
    name = :me
    m = samplesystem(;name)
    @unpack p1,x = m 
    @parameters t 
    ps = @parameters RHS_0=-0.1
    D = Differential(t)
    eqs = [
        p1 ~ RHS_0 * x,  # but p1 is not a right-hand-side item
    ]
    @test_throws ErrorException sys_ext = override_system(eqs, m; name, ps) 
end;

