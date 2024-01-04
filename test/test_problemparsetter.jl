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
popt1s = CA.ComponentVector(state = CA.ComponentVector(L = 10.1),
    par = CA.ComponentVector(k_L = 1.1, k_R = 1 / 20.1))
# use Axis for type stability, but here, check with non-typestable ps
#ps = @inferred ODEProblemParSetter(Axis(keys(u1)),Axis(keys(p1)),Axis(keys(popt)))
ps = ps1 = ODEProblemParSetter(u1, p1, popt1)

# entries with substructure
u1c = CA.ComponentVector(a = (a1 = 1, a2 = (a21 = 21, a22 = 22.0)))
p1c = CA.ComponentVector(b = (b1 = 0.1, b2 = 0.2), c = [0.01, 0.02], d = 3.0)
# as long as ComponentArrays does not support Axis-indexing, focus on top-level components rather than implementing this indexing in MTKHelpers
# note: no a1 no b, 
# a2 and c need to have correct length for updating
#poptc = CA.ComponentVector(a=(a2=1:2,), c=1:2) 
poptc = vcat(u1c[CA.KeepIndex(:a)], p1c[(:b, :c)])
poptcs = CA.ComponentVector(state = u1c[CA.KeepIndex(:a)], par = p1c[(:b, :c)])
psc = pset = ODEProblemParSetter(u1c, p1c, poptc)
#u0 = u1c; p=p1c; popt=poptc

# test states and parameters and CA.ComponentVector{SA.SVector}
u1s = label_state(psc, SA.SVector{3}(CA.getdata(u1c)))
p1s = label_par(psc, SA.SVector{5}(CA.getdata(p1c))) # convert to CA.ComponentVector{SA.SVector}

@testset "_get_axis of ComponentVectors and Strings" begin
    popt_strings_tup = ("L", "k_L", "k_R")
    ps = ODEProblemParSetter(u1, p1, popt_strings_tup)
    tmp_state = CA.ComponentVector(L = 1)
    tmp_par = CA.ComponentVector(k_L = 2, k_R = 3)
    tmp = CA.ComponentVector(state = tmp_state, par = tmp_par)
    @test axis_paropt(ps) == MTKHelpers._get_axis(tmp)
end;

@testset "axis_length FlatAxis" begin
    # empty substructure gives a FlatAxis
    u1_flat = CA.ComponentVector{Float64}()
    popt1_no_u = CA.ComponentVector(k_L = 1.1, k_R = 1 / 20.1)
    ps = ODEProblemParSetter(u1_flat, p1, popt1_no_u)
    @test count_state(ps) == 0
end;

@testset "_ax_symbols" begin
    cv = CA.ComponentVector(a = (a1 = 100, a2 = (a21 = 210, a22 = 220)),
        c = (c1 = reshape(1:4, (2, 2)),))
    ax = first(CA.getaxes(cv))
    #tp = @inferred MTKHelpers._ax_symbols_tuple(ax; prefix = "_") # lastindex not inferable
    tp = MTKHelpers._ax_symbols_tuple(ax; prefix = "_")
    @test tp == (:a_a1,
        :a_a2_a21,
        :a_a2_a22,
        Symbol("c_c1[1,1]"),
        Symbol("c_c1[2,1]"),
        Symbol("c_c1[1,2]"),
        Symbol("c_c1[2,2]"))
    v = @inferred MTKHelpers._ax_symbols_vector(ax, prefix = "_")
    @test v == collect(tp)
end

@testset "warning on some optimization parameters not found in u0 nor p" begin
    state_syms = keys(u1)
    par_syms = keys(p1)
    popt_syms = (:L, :k_L, :M1, :M2)
    @test_throws ErrorException ODEProblemParSetter(state_syms, par_syms, popt_syms)
    psw = @test_logs (:warn, r"M1.+M2") ODEProblemParSetter(CA.Axis(state_syms),
        CA.Axis(par_syms), CA.Axis(popt_syms))
    #get_paropt(psw, u1, p1) # error, because Missing not allowed
    #test if setting parameters does work
end;

@testset "warning on duplicate symbols" begin
    state_syms = keys(p1)
    par_syms = keys(p1)
    popt_syms = (:k_L,)
    @test_throws ErrorException ODEProblemParSetter(state_syms, par_syms, popt_syms)
    psw = @test_logs (:warn, r"k_L") ODEProblemParSetter(CA.Axis(state_syms),
        CA.Axis(par_syms), CA.Axis(popt_syms))
    p2 = p1 .* 2

    symbols_paropt(psw)
    res = get_paropt_labeled(psw, p2, p1)
    res = get_paropt(psw, p2, p1)
    @test res == [p2[1]] # picked the state (u0) value
    #
    # specify a different underlying type
    # res = get_paropt(psw, p2, p1, MTKHelpers.ODEMVectorCreator())
    # @test_broken res isa MVector
end;

@testset "get_paropt_labeled with StaticVector" begin
    # extract names and values from SA.SVector
    pset = ODEProblemParSetter(u1s, p1s, poptc)
    res_vec = get_paropt(pset, u1s, p1s)
    @test res_vec isa AbstractVector{Float64}
end;

# @testset "MethodError if missing symbols are not allowed" begin
#     @test_throws MethodError ps1 = ODEProblemParSetter((:x, :RHS), (:τ,), (:RHS, :τ, :M), Val(false))
# end;

@testset "access keys and counts" begin
    # axis is not inferred (AbstractAxis)
    # @test keys((@inferred axis_state(ps))) == keys(u1)
    # @test keys((@inferred axis_par(ps))) == keys(p1)
    # @test keys((@inferred axis_paropt(ps))) == (:state,:par)
    # also not with function barrier
    # tmpf = (ps) -> begin
    #     @test keys((@inferred axis_state(ps))) == keys(u1)
    #     @test keys((@inferred axis_par(ps))) == keys(p1)
    #     @test keys((@inferred axis_paropt(ps))) == (:state,:par)
    # end
    # tmpf(ps)
    @test keys((axis_state(ps))) == keys(u1)
    @test keys((axis_par(ps))) == keys(p1)
    @test keys((axis_paropt(ps))) == (:state, :par)
    #
    @test (@inferred count_state(ps)) == length(u1)
    @test (@inferred count_par(ps)) == length(p1)
    @test (@inferred count_paropt(ps)) == length(popt1)
end;

@testset "access symbols" begin
    # u1c = CA.ComponentVector(a=(a1=1,a2=(a21=21, a22=22.0)))
    # p1c = CA.ComponentVector(b=(b1=0.1, b2=0.2), c=[0.01, 0.02], d=3.0)
    # as long as ComponentArrays does not support Axis-indexing, focus on top-level components rather than implementing this indexing in MTKHelpers
    # note: no a1 no b, 
    # a2 and c need to have correct length for updating
    #poptc = CA.ComponentVector(a=(a2=1:2,), c=1:2) 
    poptc = vcat(u1c[CA.KeepIndex(:a)], p1c[(:b, :c)])
    @test (symbols_state(psc)) == (:a₊a1, :a₊a2₊a21, :a₊a2₊a22)
    @test (symbols_par(psc)) == (:b₊b1, :b₊b2, Symbol("c[1]"), Symbol("c[2]"), :d)
    # symbols paropt should neglect classification of symbols, which can be 
    # obtained by keys(axis_paropt(ps))
    # @test (symbols_paropt(psc)) ==
    #       (:state₊a₊a1, :state₊a₊a2₊a21, :state₊a₊a2₊a22, :par₊b₊b1, :par₊b₊b2, Symbol("par₊c[1]"), Symbol("par₊c[2]"))
    @test (symbols_paropt(psc)) ==
          (:a₊a1, :a₊a2₊a21, :a₊a2₊a22, :b₊b1, :b₊b2, Symbol("c[1]"), Symbol("c[2]"))
end;
# ps.statemap
# ps.optinfo

function test_label_svectors(pset,
        u0,
        p,
        popt,
        ::Val{NU0},
        ::Val{NP},
        ::Val{NOPT}) where {NOPT, NU0, NP}
    #Main.@infiltrate_main
    # not inferred - compare to test_problemparsettertyped
    @test label_paropt(pset, popt) == popt
    @test (label_paropt(pset, convert(Array, popt))) == popt
    @test (label_paropt(pset, SA.SVector{NOPT}(CA.getdata(popt)))) == popt
    #@test @inferred(label_paropt(pset, SA.SVector{NOPT}(CA.getdata(popt)))) == popt
    @test (label_paropt(pset, SA.SVector{NOPT}(CA.getdata(popt))) |> CA.getdata) isa
          SA.SVector
    #
    @test label_state(pset, u0) == u0
    @test (label_state(pset, convert(Array, u0))) == u0
    ls = label_state(pset, SA.SVector{NU0}(CA.getdata(u0))) # error without CA.getdata
    @test ls == u0
    @test CA.getdata(ls) isa SA.SVector
    #
    @test label_par(pset, p) == p
    @test (label_par(pset, convert(Array, p))) == p
    psv = SA.SVector{NP}(CA.getdata(p))
    @test (label_par(pset, psv)) == p
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
    u0o, po = update_statepar(pset, CA.getdata(popt), CA.getdata(u0), CA.getdata(p))
    u0o, po = update_statepar(pset, popt, u0, p)
    #return u0o, po
    #@btime update_statepar($ps, $popt, $u1, $p1) 
    @test CA.getaxes(u0o) == CA.getaxes(u0)
    @test CA.getaxes(po) == CA.getaxes(p)
    @test typeof(CA.getdata(u0o)) == typeof(CA.getdata(u0))
    @test typeof(CA.getdata(po)) == typeof(CA.getdata(p))
    @test all(u0o .≈ u0_target)
    @test all(po .≈ p_target)
    #
    popt2n = get_paropt_labeled(pset, u0o, po)
    @test popt2n == popt
    #
    # Main.@infiltrate_main
    #popt2 = @inferred get_paropt(pset, u0o, po)  # typeof vector not defined
    popt2 = get_paropt(pset, u0o, po)
    @test all(popt2 .== popt)
    _ = @inferred eltype(get_paropt(pset, u0o, po))
    #
    popt2m = get_paropt_labeled(pset, collect(u0o), collect(po))
    @test popt2m == popt
end;

@testset "update_statepar vector unstructured" begin
    u1t = CA.ComponentVector(L = 10.1)
    pt = CA.ComponentVector(k_L = 1.1, k_R = 1 / 20.1, m = 2.0)
    pset = ps
    u0 = u1
    p = p1
    popt = popt1s
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

@testset "update ODEProblem optimized parameters" begin
    f = (u, p, t) -> p[1] * u
    u0 = (u1 = 1 / 2,)
    p = (p1 = 1.1, p2 = 2.0)
    tspan = (0.0, 1.0)
    #prob = ODEProblem(f,SA.SVector(Tuple(u0)),tspan,SA.SVector(Tuple(p))) # SA.SVector broken
    prob = ODEProblem(f, collect(u0), tspan, collect(p))
    #sol = solve(prob, Tsit5())
    #sol[end]
    pset = ODEProblemParSetter(keys(u0), keys(p), (:u1, :p2))
    popt = CA.ComponentVector(state = CA.ComponentVector(u1 = 1 / 4),
        par = CA.ComponentVector(p2 = 1.2))
    #
    # update_statepar
    prob2 = remake(prob, popt, pset)
    @test prob2.u0[1] == popt.state.u1
    @test prob2.p[1] == p.p1 # not updated
    @test prob2.p[2] == popt.par.p2
    #
    # get_paropt
    @test get_paropt(pset, prob2) == CA.getdata(popt)
    @test get_paropt_labeled(pset, prob2) == popt
    @test label_par(pset, prob.p) == CA.ComponentVector(p)
    #
    @test names(name_paropt(pset, prob2), 1) == [:u1, :p2]
end;

@testset "gradient _update_cv_top" begin
    cv = CA.ComponentVector(a = 1:2, b = 2, c = 3.0)
    popt = CA.ComponentVector(a = 11:12, b = 20.0)
    is_updated = SA.SVector(true, true, false)
    cv2 = MTKHelpers._update_cv_top(cv, popt)
    fcost = let cv = cv
        (popt) -> begin
            #local cv2 = @inferred MTKHelpers._update_cv_top(cv, popt)
            local cv2 = @inferred MTKHelpers._update_cv_top(cv, popt, is_updated)
            d = sum(cv2)
            d * d
        end
    end
    @inferred fcost(popt)
    ForwardDiff.gradient(fcost, popt)
    @test typeof(ForwardDiff.gradient(fcost, popt)) == typeof(popt)
    fb = (cv, popt) -> MTKHelpers._update_cv_top(cv, popt)
    res = @inferred fb(cv, popt)
end;

@testset "gradient with Vector" begin
    pset = psc
    u0 = u1c
    p = p1c
    popt = poptc
    cv = p1c
    fcost = let pset = pset, u0 = u0, p = p
        (popt) -> begin
            #local u0o, po = @inferred update_statepar(pset, popt, u0, p)
            # not inferred: run-time dispatch on computations based on u0o and po
            #   because popt.state -> Any after attach_axis
            local u0o, po = update_statepar(pset, popt, u0, p)
            E = eltype(typeof(u0o)) # u0o not inferred
            d = sum(get_paropt(pset, u0o, po))::E
            d * d
        end
    end
    # not inferred -> see concrete version pset = get_concrete(pset)
    #@inferred fcost(popt)  
    fcost(popt)
    #using Cthulhu
    #@descend_code_warntype fcost(popt)
    #@descend_code_warntype update_statepar(pset, popt, u0, p)
    # ForwardDiff.gradient not inferred
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
        u0o, po = update_statepar(pset, popt, u0, p)
        d = sum(get_paropt(pset, u0o, po))
        d * d
    end
    poptsv = SA.SVector{7}(popt)
    fcost(poptsv)
    @test typeof(ForwardDiff.gradient(fcost, poptsv)) == typeof(poptsv)
end;

@testset "gradient with AbstractVector" begin
    pset = psc
    u0 = u1c
    p = p1c
    popt = poptc
    cv = p1c
    fcost = (popt) -> begin
        u0o, po = update_statepar(pset, popt, u0, p)
        d = sum(get_paropt(pset, u0o, po))
        d * d
    end
    poptv = collect(popt)
    fcost(poptv)
    @test typeof(ForwardDiff.gradient(fcost, poptv)) == typeof(poptv)
end;

function test_system(ps1, popt_names, m)
    #Main.@infiltrate_main
    #@code_warntype ODEProblemParSetter(m, popt_names)
    # not inferred
    @test (keys(axis_state(ps1))) == (:x, :RHS)
    @test (keys(axis_par(ps1))) == (:τ, :p1, :p2, :i)
    @test (keys_paropt(ps1)) == popt_names
end
@testset "construct from ODESystem" begin
    @named m = samplesystem()
    popt_names = (:RHS, :τ)
    ps1 = ODEProblemParSetter(m, popt_names)
    # only type stable after function boundary
    test_system(ps1, popt_names, m)
end;

@testset "name_paropt" begin
    pset = ps
    xn = @inferred name_paropt(pset, collect(1:count_paropt(pset)))
    @test names(xn)[1] == collect(symbols_paropt(pset))
    @test xn[:k_R] == 3
    xn = @inferred name_state(pset, collect(1:count_state(pset)))
    @test names(xn)[1] == collect(symbols_state(pset))
    xn = @inferred name_par(pset, collect(1:count_par(pset)))
    @test names(xn)[1] == collect(symbols_par(pset))
    #
    frandsym = () -> begin
        syms_arr = rand([:L, :k_L, :k_R], 2)
        ntuple(i -> syms_arr[i], 2)
        rand() > 0.5 ? (:L, :k_L, :k_R) : (:L, :k_L)
    end
    #frandsym()
    # ?implement al{allow_missing_popt}=Val(true)
    # as long as Axis argument is passed, all the type is inferred
    #psr = @inferred ODEProblemParSetter(Axis(keys(u1)),Axis(keys(p1)),Axis(frandsym()))
    ftmp = (poptnames) -> ODEProblemParSetter(keys(u1), keys(p1), poptnames)
    #psr = @inferred ftmp(frandsym()) # not inferable: Axis is constructed from unknow syms
    psr = ftmp(frandsym()) # not inferable: Axis is constructed from unknow syms
    # use Parsetter either with explicit Axis or inside function barrier
    xl = label_paropt(psr, collect(1:count_paropt(psr))) # ok?
    xn = @inferred name_paropt(psr, collect(1:count_paropt(psr)))
end;
