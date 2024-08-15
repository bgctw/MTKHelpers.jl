#TestEnv.activate()
using Test
using MTKHelpers
using MTKHelpers: MTKHelpers as CP
using OrdinaryDiffEq, ModelingToolkit
using ComponentArrays: ComponentArrays as CA
using StaticArrays: StaticArrays as SA

pkgdir = dirname(dirname(pathof(MTKHelpers)))
include(joinpath(pkgdir, "test", "samplesystem.jl"))
include(joinpath(pkgdir, "test", "testset_utils.jl"))

#@named m = samplesystem()

(u1, p1, popt1s, prob_sys1) = get_sys_ex_scalar();
ps1 = pset = NullODEProblemParSetter(get_system(prob_sys1))
ps1c = get_concrete(ps1)
popt_empty = CA.ComponentVector(
    state = CA.ComponentVector{Float64}(), par = CA.ComponentVector())

(u1c, p1c, poptcs, prob_sys2) = get_sys_ex_vec();
poptc = flatten1(poptcs)
psc = pset = NullODEProblemParSetter(get_system(prob_sys2))
pscc = get_concrete(psc)

@testset "empty paropt" begin
    pset = ps1
    @test axis_paropt(pset) == CA.getaxes(popt_empty)[1]
    @test axis_paropt_scalar(pset) == CA.getaxes(popt_empty)[1]
    @test axis_paropt_flat1(pset) == CA.getaxes(flatten1(popt_empty))[1]
    pset = ps1c
    @test axis_paropt(pset) == CA.getaxes(popt_empty)[1]
    @test axis_paropt_scalar(pset) == CA.getaxes(popt_empty)[1]
    @test axis_paropt_flat1(pset) == CA.getaxes(flatten1(popt_empty))[1]
end;

@testset "label_par and label_state" begin
    pset = psc
    @test label_par(pset, p1c) == p1c
    @test label_state(pset, u1c) == u1c
    pset = pscc
    @test label_par(pset, p1c) == p1c
    @test label_state(pset, u1c) == u1c
end;

@testset "access keys and counts" begin
    ftest = (ps) -> begin
        @test keys((axis_state(ps))) == keys(u1)
        @test keys((axis_par(ps))) == keys(p1)
        @test keys((axis_paropt(ps))) == (:state, :par)
        #
        @test (@inferred count_state(ps)) == length(u1)
        @test (@inferred count_par(ps)) == length(p1)
        @test (@inferred count_paropt(ps)) == 0
    end
    ps = ps1
    ftest(ps)
    ps = ps1c
    ftest(ps)
end;

@testset "access symbols" begin
    ftest = (ps) -> begin
        poptc = vcat(u1c[CA.KeepIndex(:a)], p1c[(:b, :c)])
        sym_state = Symbol.(("(a(t))[1]", "(a(t))[2]", "(a(t))[3]"))
        @test (symbols_state(psc)) == sym_state
        sym_bc = Symbol.(("b[1]", "b[2]", "c[1]", "c[2]"))
        @test (symbols_par(psc)) == (sym_bc..., :d)
        @test symbols_paropt(psc) == NTuple{0, Symbol}()
    end
    ps = ps1
    ftest(ps)
    ps = ps1c
    ftest(ps)
end;
# ps.statemap
# ps.optinfo

@testset "get_par_labeled" begin
    ftest = (pset) -> begin
        @test get_par(pset, prob_sys2) == CA.getdata(p1c)
        @test get_par_labeled(pset, prob_sys2) == p1c
        @test get_state(pset, prob_sys2) == CA.getdata(u1c)
        @test get_state_labeled(pset, prob_sys2) == u1c
        @test get_paropt(pset, prob_sys2) == Float64[]
        @test get_paropt_labeled(pset, prob_sys2) == popt_empty
        @test label_paropt_flat1(pset, get_paropt(pset, prob_sys2)) ==
              CA.ComponentVector{Float64}()
    end
    pset = psc
    ftest(pset)
    pset = pscc
    ftest(pset)
end

@testset "remake scalars" begin
    pset = ps1
    ftest = (pset) -> begin
        probo = remake(prob_sys1, [], pset)
        @test get_par_labeled(pset, probo) == get_par_labeled(pset, prob_sys1)
        @test get_state_labeled(pset, probo) == get_state_labeled(pset, prob_sys1)
    end
    ftest(pset)
    pset = ps1c
    ftest(pset)
end;

@testset "get_state, get_par, get_paropt vector prob" begin
    pset = psc
    ftest = (pset) -> begin
        @test get_state(pset, prob_sys2) == collect(u1c)
        @test get_state_labeled(pset, prob_sys2) == u1c
        @test get_par(pset, prob_sys2) == collect(p1c)
        @test get_par_labeled(pset, prob_sys2) == p1c
        @test get_paropt(pset, prob_sys2) == CA.getdata(popt_empty)
        @test get_paropt_labeled(pset, prob_sys2) == popt_empty
    end
    ftest(pset)
    pset = pscc
    ftest(pset)
end;

@testset "remake vector prob" begin
    ftest = (pset) -> begin
        probo = remake(prob_sys2, Float64[], pset)
        @test get_par_labeled(pset, probo) == get_par_labeled(pset, prob_sys2)
        @test get_state_labeled(pset, probo) == get_state_labeled(pset, prob_sys2)
    end
    pset = psc
    ftest(pset)
    pset = pscc
    ftest(pset)
end;
