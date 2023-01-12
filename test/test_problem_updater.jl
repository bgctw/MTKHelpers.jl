f = (u,p,t) -> p[1]*u
u0 = (L=1/2,)
p = (k_L = 1.0, k_R = 2.0, k_P = 3.0)
tspan = (0.0,1.0)
#prob = ODEProblem(f,SVector(Tuple(u0)),tspan,SVector(Tuple(p))) # SVector broken
prob = ODEProblem(f,collect(u0),tspan,collect(p))

# update k_R and k_P based on K_L in prob
pset = ProblemParSetter(Axis(keys(u0)),Axis(keys(p)),Axis((:k_R,:k_P)))


@testset "KeysProblemParGetter" begin
    sk = (:k_L,:k_L)
    pu = ProblemUpdater(KeysProblemParGetter(sk), pset)
    prob2 = pu(prob)
    @test label_par(pset, prob2.p).k_R == p.k_L
    @test label_par(pset, prob2.p).k_P == p.k_L
end;

@testset "DummyParGetter" begin
    pu = ProblemUpdater(CP.DummyParGetter(), pset)
    prob2 = pu(prob)
    @test prob2.p == [1.0, 1.0, 10.0]
end;
