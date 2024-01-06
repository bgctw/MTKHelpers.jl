using Test
using MTKHelpers
using MTKHelpers: MTKHelpers as CP
using OrdinaryDiffEq, ModelingToolkit
using ComponentArrays: ComponentArrays as CA
#using StaticArrays: StaticArrays as SA

f = (u, p, t) -> p[1] * u
u0 = (L = 1 / 2,)
p = (k_L = 1.0, k_R = 2.0, k_P = 3.0)
tspan = (0.0, 1.0)
#prob = ODEProblem(f,SVector(Tuple(u0)),tspan,SVector(Tuple(p))) # SVector broken
prob = ODEProblem(f, collect(u0), tspan, collect(p))

@testset "KeysProblemParGetter" begin
    mapping = (:k_L => :k_R, :k_L => :k_P)
    pg = KeysProblemParGetter(mapping, keys(u0))
    # pset = get_concrete(ODEProblemParSetter(Axis(keys(u0)),Axis(keys(p)),Axis(keys(pg))))
    # pu = ProblemUpdater(pg, pset)
    pu = get_ode_problemupdater(pg, keys(u0), keys(p))
    prob2 = pu(prob)
    @test label_par(par_setter(pu), prob2.p).k_R == p.k_L
    @test label_par(par_setter(pu), prob2.p).k_P == p.k_L
    @test label_state(par_setter(pu), prob2.u0) == CA.ComponentVector(u0)
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

@testset "KeysProblemParGetter_arr" begin
    f = (u, p, t) -> p[1] * u
    u0 = CA.ComponentVector(L = 1 / 2)
    p = CA.ComponentVector(k_L = 1.0, k_R = [2.0, 3.0], k_P = [4.0, 5.0], k_L2 = 6.0)
    tspan = (0.0, 1.0)
    prob = ODEProblem(f, CA.getdata(u0), tspan, CA.getdata(p); tspan = (0.0, 1.0))
    #
    mapping = (:k_L => :k_L2, :k_R => :k_P)
    pu = get_ode_problemupdater(KeysProblemParGetter(mapping, keys(u0)), u0, p)
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
    @test_throws MethodError KeysProblemParGetter(sk, dk, keys(u0))
end;

@testset "DummyParGetter" begin
    pg = CP.DummyParGetter()
    pu = get_ode_problemupdater(pg, keys(u0), keys(p))
    prob2 = pu(prob)
    @test prob2.p == [1.0, 1.0, 10.0] # specific to DummyParGetter in problem_updater.jl
end;

@testset "get_concrete ProblemParUpdater used in cost function" begin
    mapping = (:k_L => :k_R, :k_L => :k_P)
    #ps = CA.ComponentVector(SVector{length(p)}(1:length(p)), Axis(keys(p)))
    pg = KeysProblemParGetter(mapping, keys(u0)) #not type-stable extractions
    pu = get_ode_problemupdater(pg, keys(u0), keys(p))
    @inferred pg(pu, prob) # type stable getter
    #pg = ODEKeysProblemParGetter(keys(u0), keys(p), mapping) 
    puc1 = get_concrete(pu)
    prob2 = pu(prob) # not inferred
    prob3 = @inferred puc1(prob)
    # using Cthulhu
    # @descend_code_warntype puc1(prob)
    @test prob3.p == prob2.p
    @test prob3.u0 == prob2.u0
    #
    get_fopt = (pu) -> begin
        # get a concrete-type version of the ProblemParSetter and pass it 
        # through a function barrier to a closure (function within let)
        puc = get_concrete(pu)
        get_fopt_inner = (puc) -> begin
            let puc = puc
                (prob) -> begin
                    prob_upd = @inferred puc(prob)
                end # function
            end # let
        end # get_fopt_inner  
        get_fopt_inner(puc)
    end # get_ftopt
    fopt = get_fopt(pu)
    res = @inferred fopt(prob)
    @test res.p == prob2.p
end;

@testset "ProblemParUpdater from ODESystem" begin
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
    @test label_par(pset, prob.p).m₊τ == 2.0
    prob2 = pu(prob)
    @test label_par(pset, prob2.p).m₊τ == 1.0
end;
