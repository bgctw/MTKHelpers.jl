using Test
using MTKHelpers
using MTKHelpers: MTKHelpers as CP
using OrdinaryDiffEq, ModelingToolkit
using ComponentArrays: ComponentArrays as CA
using SymbolicIndexingInterface: SymbolicIndexingInterface as SII
#using StaticArrays: StaticArrays as SA

pkg_dir = dirname(dirname(pathof(MTKHelpers)))
include(joinpath(pkg_dir, "test", "testset_utils.jl"))
include(joinpath(pkg_dir, "test", "samplesystem.jl"))


@testset "KeysProblemParGetter" begin
    (u0, p, popt, prob) = get_sys_ex_scalar();
    mapping = (:k_L => :k_R, :k_L => :m)
    pg = KeysProblemParGetter(mapping, keys(u0))
    # pset = get_concrete(ODEProblemParSetter(Axis(keys(u0)),Axis(keys(p)),Axis(keys(pg))))
    # pu = ProblemUpdater(pg, pset)
    pu = get_ode_problemupdater(pg, get_system(prob))
    prob2 = pu(prob)
    pnew = get_par_labeled(par_setter(pu), prob2)
    @test pnew.k_R == pnew.k_L
    @test pnew.m == pnew.k_L
    @test get_state_labeled(par_setter(pu), prob2) == get_state_labeled(par_setter(pu), prob)
end;

@testset "NullProblemUpdater" begin
    pu = NullProblemUpdater()
    @test get_concrete(pu) === pu
end;

@testset "no_concrete" begin
    # test warning if there is no concrete version of a problemupdater
    struct DummyProblemUpdater <: AbstractProblemUpdater end
    (pu::DummyProblemUpdater)(prob) = prob
    MTKHelpers.isconcrete(::DummyProblemUpdater) = false
    pu = DummyProblemUpdater()
    puc = @test_logs (:warn, r"DummyProblemUpdater") get_concrete(pu)
    @test puc === pu
end;

@testset_skip "KeysProblemParGetter_arr" begin
    (u0, p, poptcs, prob) = get_sys_ex_vec();
    poptc = flatten1(poptcs)
    # f = (u, p, t) -> p[1] * u
    # u0 = CA.ComponentVector(L = 1 / 2)
    # p = CA.ComponentVector(k_L = 1.0, k_R = [2.0, 3.0], k_P = [4.0, 5.0], k_L2 = 6.0)
    # tspan = (0.0, 1.0)
    # prob = ODEProblem(f, CA.getdata(u0), tspan, CA.getdata(p); tspan = (0.0, 1.0))
    #
    mapping = (:b => :c,)
    pu = get_ode_problemupdater(KeysProblemParGetter(mapping, keys(u0)), get_system(prob))
    #axis_par(par_setter(pu))
    prob2 = pu(prob)
    pset = par_setter(pu)
    p2 = get_par_labeled(pset, prob2)
    u02 = get_state_labeled(pset, prob2)
    @test p2.c == p.b
    dest = Tuple(last(p) for p in mapping)
    u0_non_opt = setdiff(keys_state(pset), dest)
    @test u02[u0_non_opt] == u0[u0_non_opt]
    p_non_opt = setdiff(keys_par(pset), dest)
    @test p2[p_non_opt] == p[p_non_opt]
end;

@testset "KeysProblemParGetter error on duplicated destinations" begin
    (u0, p, popt, prob) = get_sys_ex_scalar();
    mapping = (:k_L => :m, :k_R => :m)
    @test_throws ErrorException KeysProblemParGetter(mapping, keys(u0))
end;

@testset "DummyParGetter" begin
    (u0, p, popt, prob) = get_sys_ex_scalar();
    pg = CP.DummyParGetter()
    pu = get_ode_problemupdater(pg, get_system(prob))
    prob2 = pu(prob)
    p2 = get_par_labeled(par_setter(pu), prob2)
    @test p2[:k_L] == p[:k_L] == p2[:k_R] == p2[:m]/10 # specific to DummyParGetter in problem_updater.jl
end;

@testset "get_concrete ProblemParUpdater used in cost function" begin
    (u0, p, popt, prob) = get_sys_ex_scalar();
    mapping = (:k_L => :k_R, :k_L => :m)
    #ps = CA.ComponentVector(SVector{length(p)}(1:length(p)), Axis(keys(p)))
    pg = KeysProblemParGetter(mapping, keys(u0)) #not type-stable extractions
    pu = get_ode_problemupdater(pg, get_system(prob))
    #@inferred pg(pu, prob) # type stable getter
    #pg = ODEKeysProblemParGetter(keys(u0), keys(p), mapping) 
    puc1 = get_concrete(pu)
    prob2 = pu(prob) # not inferred
    #prob3 = @inferred puc1(prob)
    @test_broken "@inferred puc1(prob)" == "remake not inferred"
    prob3 = puc1(prob)
    # using Cthulhu
    # @descend_code_warntype puc1(prob)
    @test prob3.p == prob2.p
    @test prob3.u0 == prob2.u0
    #
    m_getter = SII.getp(prob, :m)
    @inferred m_getter(prob)
    @inferred m_getter(prob2)
    get_fopt = (pu) -> begin
        # get a concrete-type version of the ProblemParSetter and pass it 
        # through a function barrier to a closure (function within let)
        puc = get_concrete(pu)
        get_fopt_inner = (puc) -> begin
            let puc = puc
                (prob) -> begin
                    #prob_upd = @inferred puc(prob)
                    prob_upd = puc(prob)
                    m_getter(prob_upd)
                end # function
            end # let
        end # get_fopt_inner  
        get_fopt_inner(puc)
    end # get_ftopt
    fopt = get_fopt(pu)
    res = fopt(prob)
    #res = @inferred fopt(prob)
    @test res == m_getter(prob2)
end;

@testset_skip "ProblemParUpdater from ODESystem" begin
    @named m = MTKHelpers.samplesystem_vec()
    @named sys = embed_system(m)
    pset = ODEProblemParSetter(sys, (:m₊p,))
    mapping = (:m₊i => :m₊τ,)
    pg = KeysProblemParGetter(mapping, keys_state(pset))
    pu = get_ode_problemupdater(pg, sys)
    #
    st = Symbolics.scalarize(m.x .=> [1.0, 2.0])
    p_new = vcat(m.i => 1.0, m.τ => 2.0, Symbolics.scalarize(m.p .=> [2.1, 2.2, 2.3]))
    prob = ODEProblem(sys, st, (0.0, 10.0), p_new)
    @test get_par_labeled(pset, prob).m₊τ == 2.0
    prob2 = pu(prob)
    @test get_par_labeled(pset, prob2).m₊τ == 1.0
end;
