@testset "getlast" begin
    @named m = samplesystem()
    @named sys = embed_system(m)
    prob = ODEProblem(sys, [m.x => 0.0], (0.0, 10.0), [m.τ => 3.0])
    sol = solve(prob, Tsit5());
    res = @inferred getlast(sol, m.x, m.RHS)
    @test res == NamedArrays.NamedArray([sol[m.x, end], sol[m.RHS, end]], ([m.x, m.RHS],))
    res2 = @inferred getlast(sol, res) # providing NamedArray to getlast
    @test res2 == res
end;
