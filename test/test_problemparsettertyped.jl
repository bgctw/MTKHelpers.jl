#using Infiltrator

# states and parameters are single entries
u1 = ComponentVector(L = 10.0)
p1 = ComponentVector(k_L = 1.0, k_R = 1 / 20, m = 2.0)
popt1 = ComponentVector(L = 10.1, k_L = 1.1, k_R = 1 / 20.1)
popt1s = ComponentVector(state=ComponentVector(L = 10.1), 
    par=ComponentVector(k_L = 1.1, k_R = 1 / 20.1))
# use Axis for type stability, but here, check with non-typestable ps
#ps = @inferred ODEProblemParSetterTyped(Axis(keys(u1)),Axis(keys(p1)),Axis(keys(popt)))
ps = ps1 = ODEProblemParSetterTyped(u1, p1, popt1)

# entries with substructure
u1c = ComponentVector(a = (a1 = 1, a2 = (a21 = 21, a22 = 22.0)))
p1c = ComponentVector(b = (b1 = 0.1, b2 = 0.2), c = [0.01, 0.02], d = 3.0)
# as long as ComponentArrays does not support Axis-indexing, focus on top-level components rather than implementing this indexing in MTKHelpers
# note: no a1 no b, 
# a2 and c need to have correct length for updating
#poptc = ComponentVector(a=(a2=1:2,), c=1:2) 
poptc = vcat(u1c[KeepIndex(:a)], p1c[(:b, :c)])
poptcs = ComponentVector(state=u1c[KeepIndex(:a)], par=p1c[(:b, :c)])
psc = pset = ODEProblemParSetterTyped(u1c, p1c, poptc)
#u0 = u1c; p=p1c; popt=poptc

# test states and parameters and ComponentVector{SVector}
u1s = label_state(psc, SVector{3}(getdata(u1c)))
p1s = label_par(psc, SVector{5}(getdata(p1c))) # convert to ComponentVector{SVector}

@testset "access keys and counts" begin
    @test (@inferred keys(axis_state(ps))) == keys(u1)
    @test (@inferred keys(axis_par(ps))) == keys(p1)
    @test (@inferred keys(axis_paropt(ps))) == (:state,:par)
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
    @test @inferred(label_paropt(pset, SVector{NOPT}(getdata(popt)))) == popt
    @test (label_paropt(pset, SVector{NOPT}(getdata(popt))) |> getdata) isa SVector
    #
    @test label_state(pset, u0) == u0
    @test @inferred(label_state(pset, convert(Array, u0))) == u0
    ls = @inferred label_state(pset, SVector{NU0}(getdata(u0))) # error without getdata
    @test ls == u0
    @test getdata(ls) isa SVector
    #
    @test label_par(pset, p) == p
    @test @inferred(label_par(pset, convert(Array, p))) == p
    psv = SVector{NP}(getdata(p))
    @test @inferred(label_par(pset, psv)) == p
    #@btime label_par($ps, $psv) # 3 allocations? creating views for subectors
    lp = label_par(pset, psv)
    @test getdata(lp) === psv
    @test getaxes(lp) === getaxes(p)
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
    u0o, po = @inferred update_statepar(pset, getdata(popt), getdata(u0), getdata(p))
    u0o, po = @inferred update_statepar(pset, popt, u0, p)
    #return u0o, po
    #@btime update_statepar($ps, $popt, $u1, $p1) 
    @test getaxes(u0o) == getaxes(u0)
    @test getaxes(po) == getaxes(p)
    @test typeof(getdata(u0o)) == typeof(getdata(u0))
    @test typeof(getdata(po)) == typeof(getdata(p))
    @test all(u0o .≈ u0_target)
    @test all(po .≈ p_target)
    #
    #using Cthulhu
    #@descend_code_warntype get_paropt_labeled(ps, u0o, po)
    #inferred only works with CA.getdata 
    popt2n = @inferred get_paropt_labeled(pset, u0o, po)
    @test popt2n == popt
    #
    popt2 = @inferred get_paropt(pset, u0o, po)
    @test all(popt2 .== popt)
    #
    popt2m = @inferred get_paropt_labeled(pset, collect(u0o), collect(po))
    @test popt2m == popt
end;

@testset "update_statepar vector unstructured" begin
    u1t = ComponentVector(L = 10.1)
    pt = ComponentVector(k_L = 1.1, k_R = 1 / 20.1, m = 2.0)
    pset = ps
    u0 = u1
    p = p1
    popt = popt1s
    # inferred although pset is global
    _, _ = @inferred update_statepar(pset, getdata(popt), getdata(u0), getdata(p))
    _ = @inferred get_paropt_labeled(pset, collect(u0), collect(p))
    #@code_warntype get_paropt_labeled(pset, collect(u0), collect(p))
    #@descend_code_warntype get_paropt_labeled(pset, collect(u0), collect(p))
    test_update_statepar_and_get_paropt(pset, u1, p1, popt, u1t, pt)
end;
@testset "update_statepar vector structured" begin
    #u1t = ComponentVector(a = (a1 = 1, a2 = (a21 = 1, a22 = 2.0))) # only subcomponent
    u1t = ComponentVector(a = (a1 = 1, a2 = (a21 = 21, a22 = 22.0)))
    pt = ComponentVector(b = (b1 = 0.1, b2 = 0.2), c = [0.01, 0.02], d = 3.0)
    pset = psc
    u0 = u1c
    p = p1c
    popt = poptcs
    cv = p1c
    @inferred get_paropt_labeled(pset, u0, p)
    test_update_statepar_and_get_paropt(psc, u1c, p1c, poptcs, u1t, pt)
end;
@testset "update_statepar Svector structured" begin
    u1t = ComponentVector(a = (a1 = 1, a2 = (a21 = 21, a22 = 22.0)))
    pt = ComponentVector(b = (b1 = 0.1, b2 = 0.2), c = [0.01, 0.02], d = 3.0)
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
#     s = ComponentVector(k_L = 1.1, k_R = 1/20.1)
#     p1a_l = label_par(ps, p1a)
#     tmp = _update_cv(p1a_l, s)
#     u1t = ComponentVector(L = 10.1)
#     pt = ComponentVector(k_L = 1.1, k_R = 1/20.1, m = 2.0)
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
            d = sum(get_paropt(pset, u0o, po))
            d * d
        end
    end
    @inferred fcost(popt)
    # ForwardDiff.gradient is not inferred
    # ftmp = (popt) -> @inferred ForwardDiff.gradient(fcost, popt)
    # res = ftmp(popt)
    res = ForwardDiff.gradient(fcost, popt)
    @test typeof(res) == typeof(popt)
end;

@testset "gradient with SVector" begin
    pset = psc
    u0 = u1c
    p = p1c
    popt = poptc
    cv = p1c
    fcost = (popt) -> begin
        u0o, po = @inferred update_statepar(pset, popt, u0, p)
        d = sum(get_paropt(pset, u0o, po))
        d * d
    end
    poptsv = SVector{7}(popt)
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
        d = sum(get_paropt(pset, u0o, po))
        d * d
    end
    poptv = collect(popt)
    @inferred fcost(poptv)
    @test typeof(ForwardDiff.gradient(fcost, poptv)) == typeof(poptv)
end;

function test_system(ps1, popt_names, m)
    #Msin.@infiltrate_main
    #@code_warntype ODEProblemParSetterTyped(m, Axis(popt_names))
    #@descend_code_warntype ODEProblemParSetterTyped(m, Axis(popt_names))
    @test @inferred(keys(axis_state(ps1))) == (:x, :RHS)
    @test @inferred(keys(axis_par(ps1))) == (:τ, :p1, :p2, :i)
    @test @inferred(keys_paropt(ps1)) == popt_names
end
@testset "construct from ODESystem" begin
    @named m = samplesystem()
    popt_names = (:RHS, :τ)
    ps1 = ODEProblemParSetterTyped(m, Axis(popt_names))
    # only type stable after function boundary
    test_system(ps1, popt_names, m)
    #
    # @test_broken "stripping names of embedded system" == "currently not supported"
    # tmpf = () -> begin
    #     @named em = embed_system(m, simplify = false)
    #     #ps1e = ODEProblemParSetterTyped(em, popt_names; strip = true)
    #     test_system(ps1e, popt_names, em)
    # end
end;

