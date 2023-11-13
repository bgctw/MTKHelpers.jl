f = (u, p, t) -> p[1] * u
u0 = (L = 1 / 2,)
p = (k_L = 1.0, k_R = 2.0, k_P = 3.0)
tspan = (0.0, 1.0)
#prob = ODEProblem(f,SVector(Tuple(u0)),tspan,SVector(Tuple(p))) # SVector broken
prob = ODEProblem(f, collect(u0), tspan, collect(p))

@testset "KeysProblemParGetter" begin
    sk = (:k_L, :k_L)
    dk = (:k_R, :k_P)
    pg = KeysProblemParGetter(sk, dk)
    # pset = ODEProblemParSetter(Axis(keys(u0)),Axis(keys(p)),Axis(keys(pg)))
    # pu = ProblemUpdater(pg, pset)
    pu = get_ode_problemupdater(pg, keys(u0), keys(p))
    prob2 = pu(prob)
    @test label_par(par_setter(pu), prob2.p).k_R == p.k_L
    @test label_par(par_setter(pu), prob2.p).k_P == p.k_L
    @test label_state(par_setter(pu), prob2.u0) == ComponentVector(u0)
end;

@testset "KeysProblemParGetter error on length(source) != length(dest)" begin
    sk = (:k_L, :k_L)
    dk = (:k_R, :k_P, :k_LP)
    @test_throws MethodError KeysProblemParGetter(sk, dk)
end;

@testset "DummyParGetter" begin
    pg = CP.DummyParGetter()
    pu = get_ode_problemupdater(pg, keys(u0), keys(p))
    prob2 = pu(prob)
    @test prob2.p == [1.0, 1.0, 10.0] # specific to DummyParGetter in problem_updater.jl
end;
