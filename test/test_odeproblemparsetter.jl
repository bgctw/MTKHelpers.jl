#TestEnv.activate()
using Test
using MTKHelpers
using MTKHelpers: MTKHelpers as CP
using OrdinaryDiffEq, ModelingToolkit
using ComponentArrays: ComponentArrays as CA
using StaticArrays: StaticArrays as SA

using ForwardDiff: ForwardDiff

pkgdir = dirname(dirname(pathof(MTKHelpers)))
include(joinpath(pkgdir,"test","samplesystem.jl"))
include(joinpath(pkgdir,"test","testset_utils.jl"))

@named m = samplesystem()

(u1, p1, popt1s, prob_sys1) = get_sys_ex_scalar();
popt1 = flatten1(popt1s)
ps1ps = ps1 = ODEProblemParSetter(get_system(prob_sys1), popt1)
get_paropt_labeled(ps1, prob_sys1)

(u1c, p1c, poptcs, prob_sys2) = get_sys_ex_vec();
poptc = flatten1(poptcs)
psc = pset = ODEProblemParSetter(get_system(prob_sys2), poptc)
get_paropt_labeled(pset, prob_sys2)
label_par(pset,p1c)
label_state(pset, u1c)

# test states and parameters and CA.ComponentVector{SA.SVector}
u1s = label_state(psc, SA.SVector{3}(CA.getdata(u1c)))
p1s = label_par(psc, SA.SVector{5}(CA.getdata(p1c))) # convert to CA.ComponentVector{SA.SVector}

@testset "label_par with Parameterobject" begin
    @test_throws "get_par_labeled" label_par(ps1, prob_sys1.p)
end;

@testset "_get_axis of ComponentVectors and Strings" begin
    popt_strings_tup = ("L", "k_L", "k_R")
    ps = ODEProblemParSetter(u1, p1, popt_strings_tup, get_system(prob_sys1))
    tmp_state = CA.ComponentVector(L = 1)
    tmp_par = CA.ComponentVector(k_L = 2, k_R = 3)
    tmp = CA.ComponentVector(state = tmp_state, par = tmp_par)
    @test axis_paropt(ps) == MTKHelpers._get_axis(tmp)
end;

@testset "axis_length FlatAxis" begin
    # empty substructure gives a FlatAxis
    p1_flat = CA.ComponentVector{Float64}()
    sts = @variables L(t)
    ps = Num[]
    eq = [D(L) ~ 0, ]
    sys_flat = ODESystem(eq, t, sts, vcat(ps...); name=:sys_flat)
    popt1_flat = CA.ComponentVector(L = 10.1,)
    ps = ODEProblemParSetter(u1, p1_flat, popt1_flat, sys_flat)
    @test count_par(ps) == 0
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
    sys1 = get_system(prob_sys1)
    @test_throws ErrorException ODEProblemParSetter(state_syms, par_syms, popt_syms, sys1)
    psw = @test_logs (:warn, r"M1.+M2") (:warn,) ODEProblemParSetter(CA.Axis(state_syms),
        CA.Axis(par_syms), CA.Axis(popt_syms), sys1)
    #get_paropt(psw, u1, p1) # error, because Missing not allowed
    #test if setting parameters does work
end;


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
    ps = ps1
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
    ps = ps1
    poptc = vcat(u1c[CA.KeepIndex(:a)], p1c[(:b, :c)])
    #@test (symbols_state(psc)) == (:a₊a1, :a₊a2₊a21, :a₊a2₊a22)
    #@test (symbols_state(psc)) == (:a₊a1, :a₊a2, :a₊a3)
    sym_state = Symbol.(("(a(t))[1]", "(a(t))[2]", "(a(t))[3]"))
    @test (symbols_state(psc)) == sym_state
    #@test (symbols_state(psc)) == (Symbol("getindex(a(t), 1)"), Symbol("getindex(a(t), 2)"), Symbol("getindex(a(t), 3)"))
    #@test (symbols_par(psc)) == (:b₊b1, :b₊b2, Symbol("c[1]"), Symbol("c[2]"), :d)
    sym_bc = Symbol.(("b[1]", "b[2]","c[1]","c[2]"))
    @test (symbols_par(psc)) == (sym_bc..., :d)
    # symbols paropt should neglect classification of symbols, which can be 
    # obtained by keys(axis_paropt(ps))
    # @test (symbols_paropt(psc)) ==
    #       (:state₊a₊a1, :state₊a₊a2₊a21, :state₊a₊a2₊a22, :par₊b₊b1, :par₊b₊b2, Symbol("par₊c[1]"), Symbol("par₊c[2]"))
    @test (symbols_paropt(psc)) == (sym_state..., sym_bc...)
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
    @test label_paropt_flat1(pset, popt) == flatten1(popt)
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
    test_label_svectors(ps1, u1, p1, popt1s, Val(1), Val(3), Val(3))
end;
@testset "label Vectors structured" begin
    test_label_svectors(psc, u1c, p1c, poptcs, Val(3), Val(5), Val(7))
end;
@testset "label SVectors structured" begin
    test_label_svectors(psc, u1s, p1s, poptcs, Val(3), Val(5), Val(7))
end;

# only allowed for system, not for parameters any more
# @testset "update ODEProblem optimized parameters" begin
#     f = (u, p, t) -> p[1] * u
#     u0 = (u1 = 1 / 2,)
#     p = (p1 = 1.1, p2 = 2.0)
#     tspan = (0.0, 1.0)
#     #prob = ODEProblem(f,SA.SVector(Tuple(u0)),tspan,SA.SVector(Tuple(p))) # SA.SVector broken
#     prob = ODEProblem(f, collect(u0), tspan, collect(p))
#     #sol = solve(prob, Tsit5())
#     #sol[end]
#     pset = ODEProblemParSetter(keys(u0), keys(p), (:u1, :p2))
#     popt = CA.ComponentVector(state = CA.ComponentVector(u1 = 1 / 4),
#         par = CA.ComponentVector(p2 = 1.2))
#     #
#     # update_statepar
#     prob2 = remake(prob, popt, pset)
#     @test prob2.u0[1] == popt.state.u1
#     @test prob2.p[1] == p.p1 # not updated
#     @test prob2.p[2] == popt.par.p2
#     #
#     # get_paropt
#     @test get_paropt(pset, prob2) == CA.getdata(popt)
#     @test get_paropt_labeled(pset, prob2) == popt
#     @test label_par(pset, prob.p) == CA.ComponentVector(p)
#     #
#     @test names(name_paropt(pset, prob2), 1) == [:u1, :p2]
# end;

@testset "get_state, get_par, get_paropt scalars" begin
    @test get_state(ps1, prob_sys1) == collect(u1)
    @test get_state_labeled(ps1, prob_sys1) == u1
    @test get_par(ps1, prob_sys1) == collect(p1)
    @test get_par_labeled(ps1, prob_sys1) == p1
    popt1_orig = CA.ComponentVector(
        state=u1[keys(popt1s.state)], 
        par = p1[keys(popt1s.par)])
    @test get_paropt(ps1, prob_sys1) == collect(popt1_orig)
    @test get_paropt_labeled(ps1, prob_sys1) == popt1_orig
end;

@testset "remake scalars" begin
    probo = remake(prob_sys1, popt1s, ps1)
    @test get_paropt_labeled(ps1, probo) == popt1s
end;

@testset "get_state, get_par, get_paropt vector prob" begin
    @test get_state(psc, prob_sys2) == collect(u1c)
    @test get_state_labeled(psc, prob_sys2) == u1c
    @test get_par(psc, prob_sys2) == collect(p1c)
    @test get_par_labeled(psc, prob_sys2) == p1c
    popt1_orig = CA.ComponentVector(
        state=u1c[keys(poptcs.state)], 
        par = p1c[keys(poptcs.par)])
    @test get_paropt(psc, prob_sys2) == collect(popt1_orig)
    @test get_paropt_labeled(psc, prob_sys2) == popt1_orig
end;

@testset "remake vector prob" begin
    probo = remake(prob_sys2, poptcs, psc)
    @test get_paropt_labeled(psc, probo) == poptcs
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


@testset "gradient prob simple" begin
    prob = prob_sys1
    sys1 = get_system(prob_sys1)
    #sol = solve(prob, Tsit5())
    #sol[end]
    pset = ODEProblemParSetter(sys1, (:L, :m))
    popt = CA.ComponentVector(state = u1[(:L,)]*2, par = p1[(:m,)]*2)
    #
    # update_statepar
    #prob2 = @inferred remake(prob, popt, pset)
    prob2 = remake(prob, popt, pset)

    # setter does not work with gradient
    # setter! = SymbolicIndexingInterface.setp(sys1, [sys1.m])
    # setter!(prob2, [popt.par.m])
    # tmp = (popt) -> begin
    #     setter!(prob, popt)
    #     prob.ps[sys1.m]
    # end
    # tmp([popt.par.m])
    # res = ForwardDiff.gradient(tmp, [popt.par.m])


    u0 = u1
    p = p1
    fcost = let pset = pset, u0 = u0, p = p
        (popt) -> begin
            #local u0o, po = @inferred update_statepar(pset, popt, u0, p)
            # not inferred: run-time dispatch on computations based on u0o and po
            #   because popt.state -> Any after attach_axis
            #probo = @inferred remake(prob, popt, pset) # 
            probo = remake(prob, popt, pset)
            #local u0o, po = update_statepar(pset, popt, u0, p)
            paropt = get_paropt(pset, probo)
            #E = eltype(typeof(u0)) # u0o not inferred
            d = sum(paropt)
            d * d
        end
    end
    # not inferred -> see concrete version pset = get_concrete(pset)
    #@inferred fcost(popt)  
    fcost(popt.*3)
    #using Cthulhu
    #@descend_code_warntype fcost(popt)
    #@descend_code_warntype update_statepar(pset, popt, u0, p)
    # ForwardDiff.gradient not inferred
    res = ForwardDiff.gradient(fcost, popt)
    @test typeof(res) == typeof(popt)
    #@inferred remake(prob, popt, pset)
    #using Cthulhu
    #@descend_code_warntype remake(prob, popt, pset)
    #@descend_code_warntype update_statepar(pset, popt, u0, p)
end;

@testset "gradient with Vectors" begin
    pset = psc
    prob = prob_sys2
    #sol = solve(prob, Tsit5())
    #sol[end]
    pset = ODEProblemParSetter(get_system(prob_sys2), (:a, :b))
    popt = CA.ComponentVector(state = u1c.a.*2, par = p1c.b.*2)

    # update_statepar
    prob2 = remake(prob, popt, pset)
    u0 = u1c
    p = p1c
    cv = p1c
    fcost = let pset = pset, u0 = u0, p = p
        (popt) -> begin
            #local u0o, po = @inferred update_statepar(pset, popt, u0, p)
            # not inferred: run-time dispatch on computations based on u0o and po
            #   because popt.state -> Any after attach_axis
            probo = remake(prob, popt, pset)
            #local u0o, po = update_statepar(pset, popt, u0, p)
            paropt = get_paropt(pset, probo)
            #E = eltype(typeof(u0)) # u0o not inferred
            d = sum(paropt)
            d * d
        end
    end
    @testset "plain Vector" begin
            # not inferred -> see concrete version pset = get_concrete(pset)
        #@inferred fcost(popt)  
        fcost(popt.*3)
        #using Cthulhu
        #@descend_code_warntype fcost(popt)
        #@descend_code_warntype update_statepar(pset, popt, u0, p)
        # ForwardDiff.gradient not inferred
        res = ForwardDiff.gradient(fcost, popt)
        @test typeof(res) == typeof(popt)
    end
    @testset "gradient with SA.SVector" begin
        poptsv = SA.SVector{5}(popt)
        fcost(poptsv)
        @test typeof(ForwardDiff.gradient(fcost, poptsv)) == typeof(poptsv)
    end
end;

@testset "name_paropt" begin
    pset = ps1
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
    ftmp = (poptnames) -> ODEProblemParSetter(keys(u1), keys(p1), poptnames, get_system(prob_sys1))
    #psr = @inferred ftmp(frandsym()) # not inferable: Axis is constructed from unknow syms
    psr = ftmp(frandsym()) # not inferable: Axis is constructed from unknow syms
    # use Parsetter either with explicit Axis or inside function barrier
    xl = label_paropt(psr, collect(1:count_paropt(psr))) # ok?
    xn = @inferred name_paropt(psr, collect(1:count_paropt(psr)))
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
    popt_names = (:RHS, :τ)
    ps1 = ODEProblemParSetter(m, popt_names)
    # only type stable after function boundary
    test_system(ps1, popt_names, m)
end;

@testset "states after parameters" begin
    popt_names = (:p1, :RHS, :τ) # note parameter :p1 before state :RHS
    ps1 = @test_logs (:warn, r"state keys first") ODEProblemParSetter(m, popt_names)
    #
    # will fail:
    #test_system(ps1, popt_names, m)
end;

@testset "vcat_statesfirst" begin
    popt1 = CA.ComponentVector(p1=1)
    popt2 = CA.ComponentVector(RHS=2, τ=3)
    popt = vcat(popt1, popt2) # p1 now comes before state RHS
    popt_ordered = vcat_statesfirst(popt1, popt2; system = m)
    @test keys(popt_ordered) == (:RHS, :p1, :τ)
    @test popt_ordered[keys(popt)] == popt
end;    

@testset "warning when not flat version can be setup" begin
    function get_sys_dup()
        sts = @variables m(t)
        ps = @parameters k_L, k_R, m 
        eq = [D(m) ~ 0, ]
        sysd = ODESystem(eq, t, sts, vcat(ps...); name=:sysd)
    end
    #sysd = structural_simplify(get_sys_dup())
    sysd = complete(get_sys_dup(); split=true)
    #
    u1 = CA.ComponentVector(m = 10.0)
    p1 = CA.ComponentVector(k_L = 1.0, k_R = 1 / 20, m = 2.0)
    popt1s = CA.ComponentVector(state = u1, par = p1[(:k_L,:m)])
    ps1 = @test_logs (:warn, r"flat") ODEProblemParSetter(sysd, popt1s)
    @test (popt2 = label_paropt(ps1, CA.getdata(popt1s))) == popt1s
    # flat axis has no names, here
    @test axis_paropt_flat1(ps1) isa CA.FlatAxis
    popt2f = @test_logs (:warn, r"no flat")  label_paropt_flat1(ps1, CA.getdata(popt1s))
    @test all(popt2f .== CA.getdata(popt1s))
    #
    ps1c = get_concrete(ps1)
    popt2f = @test_logs (:warn, r"no flat")  label_paropt_flat1(ps1c, CA.getdata(popt1s))
end;



