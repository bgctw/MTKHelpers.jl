#using Infiltrator

# states and parametes are single entries
u1 = ComponentVector(L = 10.0,)
p1 = ComponentVector(k_L = 1.0, k_R = 1/20, m = 2.0)
popt1 = ComponentVector(L = 10.1, k_L = 1.1, k_R = 1/20.1)
# use Axis for type stability, but here, check with non-typestable ps
#ps = @inferred ProblemParSetter(Axis(keys(u1)),Axis(keys(p1)),Axis(keys(popt)))
ps = ps1 = ProblemParSetter(u1,p1,popt1)

# entries with substructure
u1c = ComponentVector(a=(a1=1,a2=(a21=21, a22=22.0)))
p1c = ComponentVector(b=(b1=0.1, b2=0.2), c=[0.01, 0.02], d=3.0)
# as long as ComponentArrays does not support Axis-indexing, focus on top-level compoents rather than implementing this indexing in MTKHelpers
# note: no a1 no b, 
# a2 and c need to have correct length for updating
#poptc = ComponentVector(a=(a2=1:2,), c=1:2) 
poptc = vcat(u1c[KeepIndex(:a)],p1c[(:b,:c)])
psc = pset = ProblemParSetter(u1c, p1c, poptc)
#u0 = u1c; p=p1c; popt=poptc

# test states and parameters and ComponentVector{SVector}
u1s =  label_state(psc,SVector{3}(getdata(u1c))) 
p1s =  label_par(psc,SVector{5}(getdata(p1c))) # convert to ComponentVector{SVector}


@testset "_ax_symbols" begin
    cv = ComponentVector(a=(a1=100,a2=(a21=210, a22=220)), c = (c1=reshape(1:4,(2,2)),))
    ax = first(getaxes(cv))
    #tp = @inferred MTKHelpers._ax_symbols_tuple(ax; prefix = "_") # lastindex not inferrable
    tp = MTKHelpers._ax_symbols_tuple(ax; prefix = "_") 
    @test tp == (
        :a_a1, :a_a2_a21, :a_a2_a22, 
        Symbol("c_c1[1,1]"), Symbol("c_c1[2,1]"), Symbol("c_c1[1,2]"), Symbol("c_c1[2,2]"))
    v = @inferred MTKHelpers._ax_symbols_vector(ax, prefix = "_") 
    @test v == collect(tp)
end

@testset "warning on missing symbols" begin
    state_syms = keys(u1)
    par_syms = keys(p1)
    popt_syms = (:L, :k_L, :M1, :M2)
    psw = @test_logs (:warn,r"missing optimization parameters") ProblemParSetter(Axis(state_syms), Axis(par_syms), Axis(popt_syms))
    #get_paropt(psw, u1, p1) # error, because Missing not allowed
    #test if setting parameters does work
end;

@testset "warning on duplicate symbols" begin
    state_syms = keys(p1)
    par_syms = keys(p1)
    popt_syms = (:k_L,)
    psw = @test_logs (:warn,r"to be distinct") ProblemParSetter(Axis(state_syms), Axis(par_syms), Axis(popt_syms))
    p2 = p1 .* 2
    res = get_paropt(psw, p2, p1)
    @test res == [p2[1]] # picked the state (u0) value
end;

# @testset "MethodError if missing symbols are not allowed" begin
#     @test_throws MethodError ps1 = ProblemParSetter((:x, :RHS), (:τ,), (:RHS, :τ, :M), Val(false))
# end;

@testset "access keys and counts" begin
    @test (@inferred keys(axis_state(ps))) == keys(u1)
    @test (@inferred keys(axis_par(ps))) == keys(p1)
    @test (@inferred keys(axis_paropt(ps))) == keys(popt1)
    #
    @test (@inferred count_state(ps)) == length(u1)
    @test (@inferred count_par(ps)) == length(p1)
    @test (@inferred count_paropt(ps)) == length(popt1)
end;

@testset "access symbols and counts" begin
    # u1c = ComponentVector(a=(a1=1,a2=(a21=21, a22=22.0)))
    # p1c = ComponentVector(b=(b1=0.1, b2=0.2), c=[0.01, 0.02], d=3.0)
    # as long as ComponentArrays does not support Axis-indexing, focus on top-level compoents rather than implementing this indexing in MTKHelpers
    # note: no a1 no b, 
    # a2 and c need to have correct length for updating
    #poptc = ComponentVector(a=(a2=1:2,), c=1:2) 
    poptc = vcat(u1c[KeepIndex(:a)],p1c[(:b,:c)])
        @test (symbols_state(psc)) == (:a₊a1, :a₊a2₊a21, :a₊a2₊a22)
    @test (symbols_par(psc)) == (:b₊b1, :b₊b2, Symbol("c[1]"), Symbol("c[2]"), :d)
    @test (symbols_paropt(psc)) == (
        :a₊a1, :a₊a2₊a21, :a₊a2₊a22,  :b₊b1, :b₊b2, Symbol("c[1]"), Symbol("c[2]"))
end;
# ps.statemap
# ps.optinfo

function test_label_svectors(pset, u0, p, popt, ::Val{NU0}, ::Val{NP}, ::Val{NOPT}) where {NOPT,NU0,NP}
    @test label_paropt(pset, popt) == popt
    @test @inferred(label_paropt(pset, convert(Array, popt))) == popt
    @test @inferred(label_paropt(pset, SVector{NOPT}(popt))) == popt
    @test (label_paropt(pset, SVector{NOPT}(popt)) |> getdata) isa SVector
    #
    @test label_state(pset, u0) == u0
    @test @inferred(label_state(pset, convert(Array, u0))) == u0
    tmp = label_state(pset, SVector{NU0}(u0))
    ls = @inferred label_state(pset, SVector{NU0}(u0))
    @test ls == u0
    @test getdata(ls) isa SVector
    #
    @test label_par(pset, p) == p
    @test @inferred(label_par(pset, convert(Array, p))) == p
    psv = SVector{NP}(p)
    @test @inferred(label_par(pset, psv)) == p
    #@btime label_par($ps, $psv) # 3 allocations? creating views for subectors
    lp = label_par(pset, psv)
    @test getdata(lp) === psv
    @test getaxes(lp) === getaxes(p)
end

@testset "label Vectors unstructured" begin
    test_label_svectors(ps, u1, p1, popt1, Val(1), Val(3), Val(3))
end;
@testset "label Vectors structured" begin
    test_label_svectors(psc, u1c, p1c, poptc, Val(3), Val(5), Val(7))
end;
@testset "label SVectors structured" begin
    test_label_svectors(psc, u1s, p1s, poptc, Val(3), Val(5), Val(7))
end;

function test_update_statepar_and_get_paropt(pset, u0, p, popt, u0_target, p_target) 
    u0o, po = update_statepar(pset, popt, u0, p)
    #@descend_code_warntype update_statepar(pset, popt, u0, p)
    #@code_warntype update_statepar(pset, popt, u0, p)
    u0o, po = @inferred update_statepar(pset, getdata(popt), getdata(u0), getdata(p))
    u0o, po = @inferred update_statepar(pset, popt, u0, p)
    #return u0o, po
    #@btime update_statepar($ps, $popt, $u1, $p1) # zero allocations
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
    pt = ComponentVector(k_L = 1.1, k_R = 1/20.1, m = 2.0)
    pset = ps; u0 = u1; p=p1; popt = popt1
    # inferred although pset is global
    _, _ = @inferred update_statepar(pset, getdata(popt), getdata(u0), getdata(p))
    _ = @inferred get_paropt_labeled(pset, collect(u0), collect(p))
    test_update_statepar_and_get_paropt(pset, u1, p1, popt, u1t, pt)
end;
@testset "update_statepar vector structured" begin
    #u1t = ComponentVector(a = (a1 = 1, a2 = (a21 = 1, a22 = 2.0))) # only subcomponent
    u1t = ComponentVector(a = (a1 = 1, a2 = (a21 = 21, a22 = 22.0)))
    pt = ComponentVector(b = (b1 = 0.1, b2 = 0.2), c = [0.01, 0.02], d = 3.0)
    pset = psc; u0 = u1c; p=p1c; popt = poptc; cv=p1c
    @inferred get_paropt_labeled(pset, u0, p)
    test_update_statepar_and_get_paropt(psc, u1c, p1c, poptc, u1t, pt)
end;
@testset "update_statepar Svector structured" begin
    u1t = ComponentVector(a = (a1 = 1, a2 = (a21 = 21, a22 = 22.0)))
    pt = ComponentVector(b = (b1 = 0.1, b2 = 0.2), c = [0.01, 0.02], d = 3.0)
    # fix bug in ComponentArrays where KeepIndex retuns an Vector instead SVector based 
    @test_broken test_update_statepar_and_get_paropt(psc, u1s, p1s, poptc, u1t, pt)
end
# @testset "upate_statepar and get_paropt for AxisArray" begin
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

# @testset "merge: create a mofified popt AbstractVector" begin
#     popt3 = @inferred merge(popt, (k_L = 1.2,))
#     @test length(popt3) == length(popt)
#     @test popt3.k_L == 1.2
#     @test popt3.L == popt.L
#     @test popt3.k_R == popt.k_R
# end;

# @testset "merge: create a mofified popt AbstractVector on LVector" begin
#     poptlv = LArray(popt)
#     popt3 = @inferred merge(poptlv, (k_L = 1.2,))
#     @test length(popt3) == length(popt)
#     @test popt3.k_L == 1.2
#     @test popt3.L == popt.L
#     @test popt3.k_R == popt.k_R
# end;

@testset "update ODEProblem optimized parameters" begin
    f = (u,p,t) -> p[1]*u
    u0 = (u1=1/2,)
    p = (p1=1.1,p2=2.0)
    tspan = (0.0,1.0)
    #prob = ODEProblem(f,SVector(Tuple(u0)),tspan,SVector(Tuple(p))) # SVector broken
    prob = ODEProblem(f,collect(u0),tspan,collect(p))
    #sol = solve(prob)
    #sol[end]
    pset = ProblemParSetter(Axis(keys(u0)),Axis(keys(p)),Axis((:u1,:p2)))
    popt = ComponentVector(u1=1/4, p2=1.2)
    #
    # update_statepar
    prob2 = update_statepar(pset, popt, prob)
    @test prob2.u0[1] == popt.u1
    @test prob2.p[1] == p.p1 # not updated
    @test prob2.p[2] == popt.p2
    #
    # get_paropt
    @test get_paropt(pset, prob2) == getdata(popt) 
    @test get_paropt_labeled(pset, prob2) == popt 
    @test label_par(pset, prob.p) == ComponentVector(p)
end;

@testset "gradient _update_cv_top" begin
    cv = ComponentVector(a=1:2, b=2, c=3.0)
    popt = ComponentVector(a=11:12, b=20.0)
    cv2 = MTKHelpers._update_cv_top(cv, popt)
    fcost = (popt) -> begin
        cv2 = MTKHelpers._update_cv_top(cv, popt)
        d = sum(cv2)
        d*d
    end
    fcost(popt)
    ForwardDiff.gradient(fcost, popt)
    @test typeof(ForwardDiff.gradient(fcost, popt)) == typeof(popt)
    fb = (cv, popt) -> MTKHelpers._update_cv_top(cv, popt)
    res = @inferred fb(cv, popt)
end;

@testset "gradient with Vector" begin
    pset = psc; u0 = u1c; p=p1c; popt = poptc; cv=p1c
    fcost = (popt) -> begin
        u0o, po = update_statepar(pset, popt, u0, p)
        d = sum(get_paropt(pset, u0o, po))
        d*d
    end
    fcost(popt)
    res = ForwardDiff.gradient(fcost, popt)
    @test typeof(res) == typeof(popt)
end;

@testset "gradient with SVector" begin
    pset = psc; u0 = u1c; p=p1c; popt = poptc; cv=p1c
    fcost = (popt) -> begin
        u0o, po = update_statepar(pset, popt, u0, p)
        d = sum(get_paropt(pset, u0o, po))
        d*d
    end
    poptsv = SVector{7}(popt)
    fcost(poptsv)
    @test typeof(ForwardDiff.gradient(fcost, poptsv)) == typeof(poptsv)
end;

@testset "gradient with AbstractVector" begin
    pset = psc; u0 = u1c; p=p1c; popt = poptc; cv=p1c
    fcost = (popt) -> begin
        u0o, po = update_statepar(pset, popt, u0, p)
        d = sum(get_paropt(pset, u0o, po))
        d*d
    end
    poptv = collect(popt)
    fcost(poptv)
    @test typeof(ForwardDiff.gradient(fcost, poptv)) == typeof(poptv)
end;

function test_system(ps1, popt_names, m)
    #@infiltrate
    #@code_warntype ProblemParSetter(m, Axis(popt_names))
    #@descend_code_warntype ProblemParSetter(m, Axis(popt_names))
    @test @inferred(keys(axis_state(ps1))) == (:x, :RHS)
    @test @inferred(keys(axis_par(ps1))) == (:τ,)
    @test @inferred(keys(axis_paropt(ps1))) == popt_names
end
@testset "construct from ODESystem" begin
    @named m = samplesystem()
    popt_names = (:RHS, :τ)
    ps1 = ProblemParSetter(m, Axis(popt_names))
    # only type stable after function boundary
    test_system(ps1, popt_names, m)
    #
    @named em = embed_system(m, simplify=false)
    ps1e = ProblemParSetter(em, Axis(popt_names); strip=true)
    test_system(ps1e, popt_names, em)
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
        syms_arr = rand([:L,:k_L,:k_R],2)
        ntuple(i -> syms_arr[i], 2)
        rand() > 0.5 ? (:L,:k_L,:k_R) : (:L,:k_L)
    end
    #frandsym()
    # ?implement al{allow_missing_popt}=Val(true)
    # as long as Axis argument is passed, all the type is inferred
    #psr = @inferred ProblemParSetter(Axis(keys(u1)),Axis(keys(p1)),Axis(frandsym()))
    ftmp = (poptnames) -> ProblemParSetter(Axis(keys(u1)),Axis(keys(p1)),Axis(poptnames))
    #psr = @inferred ftmp(frandsym()) # not inferrable: Axis is constructed from unknow syms
    psr = ftmp(frandsym()) # not inferrable: Axis is constructed from unknow syms
    # use Parsetter either with explicit Axis or inside function barrier
    xl = @inferred label_paropt(psr, collect(1:count_paropt(psr))) # ok?
    xn = @inferred name_paropt(psr, collect(1:count_paropt(psr)))
end;




