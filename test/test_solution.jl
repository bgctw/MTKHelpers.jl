@testset "getlast" begin
    @named m = samplesystem()
    @named sys = embed_system(m)
    prob = ODEProblem(sys, [m.x => 0.0], (0.0,10.0), [m.τ => 3.0])
    sol = solve(prob)
    res = getlast(sol, m.x, m.RHS)
    @test res == NamedArray([sol[m.x,end], sol[m.RHS,end]], ([m.x, m.RHS],))   
end;

