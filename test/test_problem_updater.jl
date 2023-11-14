f = (u, p, t) -> p[1] * u
u0 = (L = 1 / 2,)
p = (k_L = 1.0, k_R = 2.0, k_P = 3.0)
tspan = (0.0, 1.0)
#prob = ODEProblem(f,SVector(Tuple(u0)),tspan,SVector(Tuple(p))) # SVector broken
prob = ODEProblem(f, collect(u0), tspan, collect(p))

@testset "KeysProblemParGetter" begin
    mapping = (:k_L => :k_R, :k_L => :k_P)
    pg = KeysProblemParGetter(mapping)
    # pset = ODEProblemParSetter(Axis(keys(u0)),Axis(keys(p)),Axis(keys(pg)))
    # pu = ProblemUpdater(pg, pset)
    pu = get_ode_problemupdater(pg, keys(u0), keys(p))
    prob2 = pu(prob)
    @test label_par(par_setter(pu), prob2.p).k_R == p.k_L
    @test label_par(par_setter(pu), prob2.p).k_P == p.k_L
    @test label_state(par_setter(pu), prob2.u0) == ComponentVector(u0)
end;

@testset "KeysProblemParGetter_arr" begin
    f = (u,p,t) -> p[1]*u
    u0 = ComponentVector(L=1/2)
    p = ComponentVector(k_L = 1.0, k_R = [2.0,3.0], k_P = [4.0,5.0], k_L2 = 6.0)
    tspan = (0., 1.)
    prob = ODEProblem(f,getdata(u0),tspan,getdata(p);tspan = (0.,1.))
    #
    mapping = (:k_L => :k_L2, :k_R => :k_P)
    pu = get_ode_problemupdater(KeysProblemParGetter(mapping), u0, p)
    #axis_par(par_setter(pu))
    prob2 = pu(prob)
    pset = par_setter(pu)
    p2 = label_par(pset, prob2.p)
    u02 = label_state(pset, prob2.u0)
    @test p2.k_P == p.k_R
    @test p2.k_L2 == p.k_L
    dest = Tuple(last(p) for p in mapping)
    u0_non_opt = setdiff(keys_state(pset), dest)
    @test u02[u0_non_opt] == u0[u0_non_opt]
    p_non_opt = setdiff(keys_par(pset), dest)
    @test p2[p_non_opt] == p[p_non_opt]
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

