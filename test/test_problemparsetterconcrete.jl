using Test
using MTKHelpers
using MTKHelpers: MTKHelpers as CP
using OrdinaryDiffEq, ModelingToolkit
using ComponentArrays: ComponentArrays as CA
using StaticArrays: StaticArrays as SA

using ForwardDiff: ForwardDiff

test_path = splitpath(pwd())[end] == "test" ? "." : "test"
#include(joinpath(test_path,"samplesystem.jl"))
include("samplesystem.jl")


# states and parameters are single entries
u1 = CA.ComponentVector(L = 10.0)
p1 = CA.ComponentVector(k_L = 1.0, k_R = 1 / 20, m = 2.0)
popt1 = CA.ComponentVector(L = 10.1, k_L = 1.1, k_R = 1 / 20.1)
popt1s = CA.ComponentVector(state = (L = 10.1,), par = (k_L = 1.1, k_R = 1 / 20.1))
# use Axis for type stability, but here, check with non-typestable ps
#ps = @inferred ODEProblemParSetterConcrete(CA.Axis(keys(u1)),CA.Axis(keys(p1)),CA.Axis(keys(popt)))
ps = ps1 = ODEProblemParSetterConcrete(u1, p1, popt1)

# entries with substructure
u1c = CA.ComponentVector(a = (a1 = 1, a2 = (a21 = 21, a22 = 22.0)))
p1c = CA.ComponentVector(b = (b1 = 0.1, b2 = 0.2), c = [0.01, 0.02], d = 3.0)
# as long as ComponentArrays does not support Axis-indexing, focus on top-level components rather than implementing this indexing in MTKHelpers
# note: no a1 no b, 
# a2 and c need to have correct length for updating
#poptc = CA.ComponentVector(a=(a2=1:2,), c=1:2) 
poptc = vcat(u1c[CA.KeepIndex(:a)], p1c[(:b, :c)])
poptcs = CA.ComponentVector(state = u1c[CA.KeepIndex(:a)], par = p1c[(:b, :c)])
psc = pset = ODEProblemParSetterConcrete(u1c, p1c, poptc)
#u0 = u1c; p=p1c; popt=poptc

# test states and parameters and CA.ComponentVector{SA.SVector}
u1s = label_state(psc, SA.SVector{3}(CA.getdata(u1c)))
p1s = label_par(psc, SA.SVector{5}(CA.getdata(p1c))) # convert to CA.ComponentVector{SA.SVector}

@testset "access keys and counts" begin
    @test (@inferred keys(axis_state(ps))) == keys(u1)
    @test (@inferred keys(axis_par(ps))) == keys(p1)
    @test (@inferred keys(axis_paropt(ps))) == (:state, :par)
end;

function test_label_svectors(pset,
        u0,
        p,
        popt,
        ::Val{NU0},
        ::Val{NP},
        ::Val{NOPT}) where {NOPT, NU0, NP}
    #Main.@infiltrate_main
    @test label_paropt(pset, popt) == popt
    @test @inferred(label_paropt(pset, convert(Array, popt))) == popt
    @test @inferred(label_paropt(pset, SA.SVector{NOPT}(CA.getdata(popt)))) == popt
    @test (label_paropt(pset, SA.SVector{NOPT}(CA.getdata(popt))) |> CA.getdata) isa SA.SVector
    #
    @test label_state(pset, u0) == u0
    @test @inferred(label_state(pset, convert(Array, u0))) == u0
    ls = @inferred label_state(pset, SA.SVector{NU0}(CA.getdata(u0))) # error without CA.getdata
    @test ls == u0
    @test CA.getdata(ls) isa SA.SVector
    #
    @test label_par(pset, p) == p
    @test @inferred(label_par(pset, convert(Array, p))) == p
    psv = SA.SVector{NP}(CA.getdata(p))
    @test @inferred(label_par(pset, psv)) == p
    #@btime label_par($ps, $psv) # 3 allocations? creating views for subectors
    lp = label_par(pset, psv)
    @test CA.getdata(lp) === psv
    @test CA.getaxes(lp) === CA.getaxes(p)
end

@testset "label Vectors unstructured" begin
    test_label_svectors(ps, u1, p1, popt1s, Val(1), Val(3), Val(3))
end;
@testset "label Vectors structured" begin
    test_label_svectors(psc, u1c, p1c, poptcs, Val(3), Val(5), Val(7))
end;
@testset "label SVectors structured" begin
    test_label_svectors(psc, u1s, p1s, poptcs, Val(3), Val(5), Val(7))
end;

function test_update_statepar_and_get_paropt(pset, u0, p, popt, u0_target, p_target)
    #0o, po = update_statepar(pset, popt, u0, p)
    #@descend_code_warntype update_statepar(pset, popt, u0, p)
    #@code_warntype update_statepar(pset, popt, u0, p)
    u0o, po = @inferred update_statepar(pset, CA.getdata(popt), CA.getdata(u0), CA.getdata(p))
    u0o, po = @inferred update_statepar(pset, popt, u0, p)
    #return u0o, po
    #@btime update_statepar($ps, $popt, $u1, $p1) 
    @test CA.getaxes(u0o) == CA.getaxes(u0)
    @test CA.getaxes(po) == CA.getaxes(p)
    @test typeof(CA.getdata(u0o)) == typeof(CA.getdata(u0))
    @test typeof(CA.getdata(po)) == typeof(CA.getdata(p))
    @test all(u0o .≈ u0_target)
    @test all(po .≈ p_target)
    #
    #using Cthulhu
    #@descend_code_warntype get_paropt_labeled(ps, u0o, po)
    #inferred only works with CA.CA.getdata 
    popt2n = get_paropt_labeled(pset, u0o, po)
    @inferred CA.getaxes(get_paropt_labeled(pset, u0o, po))
    @inferred zeros(eltype(get_paropt_labeled(pset, u0o, po)),3)
    @test popt2n == popt
    #
    popt2 = get_paropt(pset, u0o, po)
    _ = @inferred zeros(eltype(get_paropt(pset, u0o, po)),2)
    @test all(popt2 .== popt)
    #
    popt2m = get_paropt_labeled(pset, collect(u0o), collect(po))
    _ = @inferred CA.getaxes(get_paropt_labeled(pset, collect(u0o), collect(po)))
    _ = @inferred zeros(eltype(get_paropt_labeled(pset, collect(u0o), collect(po))),3)
    @test popt2m == popt
end;

@testset "update_statepar vector unstructured" begin
    u1t = CA.ComponentVector(L = 10.1)
    pt = CA.ComponentVector(k_L = 1.1, k_R = 1 / 20.1, m = 2.0)
    pset = ps
    u0 = u1
    p = p1
    popt = popt1s
    # inferred although pset is global
    _, _ = @inferred update_statepar(pset, CA.getdata(popt), CA.getdata(u0), CA.getdata(p))
    #_ = @inferred get_paropt_labeled(pset, collect(u0), collect(p))
    #array type is variable - can be Vector, SVector, ...
    _ = @inferred CA.getaxes(get_paropt_labeled(pset, collect(u0), collect(p)))
    _ = @inferred zeros(eltype(get_paropt_labeled(pset, collect(u0), collect(p))),3)
    #@code_warntype get_paropt_labeled(pset, collect(u0), collect(p))
    #@descend_code_warntype get_paropt_labeled(pset, collect(u0), collect(p))
    test_update_statepar_and_get_paropt(pset, u1, p1, popt, u1t, pt)
end;
@testset "update_statepar vector structured" begin
    #u1t = CA.ComponentVector(a = (a1 = 1, a2 = (a21 = 1, a22 = 2.0))) # only subcomponent
    u1t = CA.ComponentVector(a = (a1 = 1, a2 = (a21 = 21, a22 = 22.0)))
    pt = CA.ComponentVector(b = (b1 = 0.1, b2 = 0.2), c = [0.01, 0.02], d = 3.0)
    pset = psc
    u0 = u1c
    p = p1c
    popt = poptcs
    cv = p1c
    @inferred CA.getaxes(get_paropt_labeled(pset, u0, p))
    @inferred zeros(eltype(get_paropt_labeled(pset, u0, p)), 3)
    test_update_statepar_and_get_paropt(psc, u1c, p1c, poptcs, u1t, pt)
end;
@testset "update_statepar Svector structured" begin
    u1t = CA.ComponentVector(a = (a1 = 1, a2 = (a21 = 21, a22 = 22.0)))
    pt = CA.ComponentVector(b = (b1 = 0.1, b2 = 0.2), c = [0.01, 0.02], d = 3.0)
    test_update_statepar_and_get_paropt(psc, u1s, p1s, poptcs, u1t, pt)
    #using BenchmarkTools
    #@btime get_paropt_labeled($psc, $u1t, $pt)
end
# @testset "update_statepar and get_paropt for AxisArray" begin
#     @test_broken "AxisArray"
#     # test with AbstractVector different from Vector
#     # _update_cv returns a Vector because _
#     u1a = AxisArray(collect(u1); row = keys(u1))
#     p1a = AxisArray(collect(p1); row = keys(p1))
#     s = CA.ComponentVector(k_L = 1.1, k_R = 1/20.1)
#     p1a_l = label_par(ps, p1a)
#     tmp = _update_cv(p1a_l, s)
#     u1t = CA.ComponentVector(L = 10.1)
#     pt = CA.ComponentVector(k_L = 1.1, k_R = 1/20.1, m = 2.0)
#     test_update_statepar_and_get_paropt(ps, u1a, p1a, popt, u1t, pt)
# end;

# @testset "merge: create a modified popt AbstractVector" begin
#     popt3 = @inferred merge(popt, (k_L = 1.2,))
#     @test length(popt3) == length(popt)
#     @test popt3.k_L == 1.2
#     @test popt3.L == popt.L
#     @test popt3.k_R == popt.k_R
# end;

# @testset "merge: create a modified popt AbstractVector on LVector" begin
#     poptlv = LArray(popt)
#     popt3 = @inferred merge(poptlv, (k_L = 1.2,))
#     @test length(popt3) == length(popt)
#     @test popt3.k_L == 1.2
#     @test popt3.L == popt.L
#     @test popt3.k_R == popt.k_R
# end;

@testset "gradient with Vector" begin
    pset = psc
    u0 = u1c
    p = p1c
    popt = poptc
    cv = p1c
    fcost = let pset = pset, u0 = u0, p = p
        (popt) -> begin
            local u0o, po = @inferred update_statepar(pset, popt, u0, p)
            v = get_paropt(pset, u0o, po)
            E = eltype(typeof(v))
            d = sum(v)::E
            d * d
        end
    end
     @inferred fcost(popt)
    # using Cthulhu
    # @descend_code_warntype fcost(popt)
    # ForwardDiff.gradient is not inferred
    # ftmp = (popt) -> @inferred ForwardDiff.gradient(fcost, popt)
    # res = ftmp(popt)
    res = ForwardDiff.gradient(fcost, popt)
    @test typeof(res) == typeof(popt)
end;

@testset "gradient with SA.SVector" begin
    pset = psc
    u0 = u1c
    p = p1c
    popt = poptc
    cv = p1c
    fcost = (popt) -> begin
        u0o, po = @inferred update_statepar(pset, popt, u0, p)
        E = eltype(typeof(u0o))
        d = sum(get_paropt(pset, u0o, po))::E
        d * d
    end
    poptsv = SA.SVector{7}(popt)
    @inferred fcost(poptsv)
    @test typeof(ForwardDiff.gradient(fcost, poptsv)) == typeof(poptsv)
end;

@testset "gradient with AbstractVector" begin
    pset = psc
    u0 = u1c
    p = p1c
    popt = poptc
    cv = p1c
    fcost = (popt) -> begin
        u0o, po = @inferred update_statepar(pset, popt, u0, p)
        E = eltype(typeof(u0o))
        d = sum(get_paropt(pset, u0o, po))::E
        d * d
    end
    poptv = collect(popt)
    @inferred fcost(poptv)
    @test typeof(ForwardDiff.gradient(fcost, poptv)) == typeof(poptv)
end;

function test_system(ps1, popt_names, m)
    #Msin.@infiltrate_main
    #@code_warntype ODEProblemParSetterConcrete(m, CA.Axis(popt_names))
    #@descend_code_warntype ODEProblemParSetterConcrete(m, CA.Axis(popt_names))
    @test @inferred(keys(axis_state(ps1))) == (:x, :RHS)
    @test @inferred(keys(axis_par(ps1))) == (:τ, :p1, :p2, :i)
    @test @inferred(keys_paropt(ps1)) == popt_names
end
@testset "construct from ODESystem" begin
    @named m = samplesystem()
    popt_names = (:RHS, :τ)
    ps1 = ODEProblemParSetterConcrete(m, CA.Axis(popt_names))
    # only type stable after function boundary
    test_system(ps1, popt_names, m)
end;

@testset "get_concrete ODEProblemParSetter used in cost function" begin
    u1 = CA.ComponentVector(L = 10.0)
    p1 = CA.ComponentVector(k_L = 1.0, k_R = 1 / 20, m = 2.0)
    popt1s = CA.ComponentVector(state = (L = 10.1,), par = (k_L = 1.1, k_R = 1 / 20.1))
    ps = ps1 = ODEProblemParSetter(u1, p1, popt1s)
    get_fopt = (ps) -> begin
        # get a concrete-type version of the ProblemParSetter and pass it 
        # through a function barrier to a closure (function within let)
        psc = get_concrete(ps)
        get_fopt_inner = (psc) -> begin
            let psc = psc
                (paropt) -> begin
                    paropt_l = @inferred label_paropt(psc, paropt)
                    (; cost = paropt_l.par.k_L + paropt_l.par.k_R, paropt_l)
                end # function
            end # let
        end # get_fopt_inner  
        get_fopt_inner(psc)
    end # get_ftopt
    fopt = get_fopt(ps)
    res = @inferred fopt(CA.getdata(popt1s) .* 2)
    @test res.cost == (popt1s.par.k_L + popt1s.par.k_R) * 2
    @test CA.getaxes(res.paropt_l) == CA.getaxes(popt1s)
end;
