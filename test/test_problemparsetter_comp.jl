using ComponentArrays
using StaticArrays, LabelledArrays
import ComponentArrays as CA
using Infiltrator

# states and parametes are single entries
u1 = ComponentVector(L = 10.0,)
p1 = ComponentVector(k_L = 1.0, k_R = 1/20, m = 2.0)
popt = ComponentVector(L = 10.1, k_L = 1.1, k_R = 1/20.1)
# use Axis for type stability, but here, check with non-typestable ps
#ps = @inferred ProblemParSetterComp1(Axis(keys(u1)),Axis(keys(p1)),Axis(keys(popt)))
ps = ProblemParSetterComp1(keys(u1),keys(p1),keys(popt))

# entries with substructure
u1c = ComponentVector(a=(a1=1,a2=(a21=21, a22=22)))
p1c = ComponentVector(b=(b1=0.1, b2=0.2), c=[0.01, 0.02])
# note: no a1 no b, 
# a2 and c need to have correct length for updating
poptc = ComponentVector(a=(a2=1:2,), c=1:2) 
psc = pset = ProblemParSetterComp1(u1c, p1c, poptc)
#u0 = u1c; p=p1c; popt=poptc

# test states and parameters and ComponentVector{SVector}
u1s =  label_state(psc,SVector{3}(getdata(u1c))) 
p1s =  label_par(psc,SVector{4}(getdata(p1c))) # convert to ComponentVector{SVector}

i_todo = () -> begin
    typeof(fixpoint(parent, u1s))
    label_par(pset,p1s)
    label_par(pset,getdata(p1s))
    res = get_paropt_labeled(pset, u1s, p1s)
    getdata(res)
    #get_paropt(pset, u1c, p1c)
    res.a.a2 = res.a.a2 .* 2
    popt = res
    s = ComponentVector(a=1:3)
    tmpf(u1s, s) = _update_cv(u1s, s)
    @inferred tmpf(u1s, s) # type stable if s is passed to function 
    @inferred ComponentVector(a=1:3) # not type stable

    fixpoint(parent,cv;fmap=typeof)
    attach_axis(typeof(fixpoint(parent,cv;fmap=typeof))(res), first(getaxes(cv)))
    convert(AT, getdata(res))

    tmpf(pset, popt, u1s, p1s) = update_statepar(pset, getdata(popt), getdata(u1s), getdata(p1s))
    tmpf(pset, popt, u1s, p1s) = update_statepar(pset, popt, u1s, p1s)
    u0new, pnew = @inferred tmpf(pset, popt, u1s, p1s)

    typeof(u1s)(typeof(getdata(u1s))(getdata(u0new)), getaxes(u1s))
    @test axes(u0new) == axes(u1c)
    @test axes(pnew) == axes(p1c)
    @test u0new.a.a2 == res.a.a2
end


@test_set "_ax_symbols" begin
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
    psw = @test_logs (:warn,r"missing optimization parameters") ProblemParSetterComp1(state_syms, par_syms, popt_syms)
    #get_paropt(psw, u1, p1) # error, because Missing not allowed
    #test if setting parameters does work
end;


@testset "warning on duplicate symbols" begin
    state_syms = keys(p1)
    par_syms = keys(p1)
    popt_syms = (:k_L,)
    psw = @test_logs (:warn,r"to be distinct") ProblemParSetterComp1(state_syms, par_syms, popt_syms)
    p2 = p1 .* 2
    res = get_paropt(psw, p2, p1)
    @test res == [p2[1]] # picked the state (u0) value
end;

# @testset "MethodError if missing symbols are not allowed" begin
#     @test_throws MethodError ps1 = ProblemParSetterComp1((:x, :RHS), (:τ,), (:RHS, :τ, :M), Val(false))
# end;

@testset "access symbols and counts" begin
    @test (@inferred symbols_state(ps)) == keys(u1)
    @test (@inferred symbols_par(ps)) == keys(p1)
    @test (@inferred symbols_paropt(ps)) == keys(popt)
    #
    @test (@inferred count_state(ps)) == length(u1)
    @test (@inferred count_par(ps)) == length(p1)
    @test (@inferred count_paropt(ps)) == length(popt)
end;

@testset "access symbols and counts" begin
    # u1c = ComponentVector(a=(a1=1,a2=(a21=21, a22=22)))
    # p1c = ComponentVector(b=(b1=0.1, b2=0.2), c=[0.01, 0.02])
    # # note: no a1 no b, 
    # # a2 and c need to have correct length for updating
    # poptc = ComponentVector(a=(a2=1:2,), c=1:2) 
    @test (@inferred symbols_state(psc)) == (:a₊a1, :a₊a2₊a21, :a₊a2₊a22)
    @test (@inferred symbols_par(psc)) == (:b₊b1, :b₊b2, Symbol("c[1]"), Symbol("c[2]"))
    @test (@inferred symbols_paropt(psc)) == (
        Symbol("a₊a2[1]"), Symbol("a₊a2[2]"), Symbol("c[1]"), Symbol("c[2]"))
    #
    @test (@inferred count_state(ps)) == length(u1)
    @test (@inferred count_par(ps)) == length(p1)
    @test (@inferred count_paropt(ps)) == length(popt)
end;
# ps.statemap
# ps.optinfo

function test_label_svectors(ps, u0, p, popt, ::Val{NU0}, ::Val{NP}, ::Val{NOPT}) where {NOPT,NU0,NP}
    @test label_paropt(ps, popt) == popt
    @test @inferred(label_paropt(ps, convert(Array, popt))) == popt
    @test @inferred(label_paropt(ps, SVector{NOPT}(popt))) == popt
    @test (label_paropt(ps, SVector{NOPT}(popt)) |> getdata) isa SVector
    #
    @test label_state(ps, u0) == u0
    @test @inferred(label_state(ps, convert(Array, u0))) == u0
    @test @inferred(label_state(ps, SVector{NU0}(u0))) == u0
    @test (label_state(ps, SVector{NU0}(u0)) |> getdata) isa SVector
    #
    @test label_par(ps, p) == p
    @test @inferred(label_par(ps, convert(Array, p))) == p
    psv = SVector{NP}(p)
    @test @inferred(label_par(ps, psv)) == p
    #@btime label_par($ps, $psv) # 3 allocations? creating views for subectors
    @test (label_par(ps, psv) |> getdata) === psv
    @test (label_par(ps, psv) |> getaxes) === getaxes(p)
end

@testset "label Vectors unstructured" begin
    test_label_svectors(ps, u1, p1, popt, Val(1), Val(3), Val(3))
end;
@testset "label Vectors structured" begin
    test_label_svectors(psc, u1c, p1c, poptc, Val(3), Val(4), Val(4))
end;
@testset "label SVectors structured" begin
    test_label_svectors(psc, u1s, p1s, poptc, Val(3), Val(4), Val(4))
end;

function test_update_statepar_and_get_paropt(ps, u0, p, popt, u0_target, p_target) 
    u0o, po = @inferred update_statepar(ps, popt, u0, p)
    #return u0o, po
    #@btime update_statepar($ps, $popt, $u1, $p1) # zero allocations
    @test getaxes(u0o) == getaxes(u0)
    @test getaxes(po) == getaxes(p)
    @test typeof(getdata(u0o)) == typeof(getdata(u0))
    @test typeof(getdata(po)) == typeof(getdata(p))
    @test all(u0o .≈ u0_target)
    @test all(po .≈ p_target)
    #
    #@descend_code_warntype get_paropt(ps, u0o, po)
    #inferred only works with CA.getdata 
    popt2 = @inferred get_paropt(ps, u0o, po)
    @test all(popt2 .== popt)
    #
    popt2n = @inferred get_paropt_labeled(ps, u0o, po)
    @test popt2n == popt
    #
    popt2m = @inferred get_paropt_labeled(ps, collect(u0o), collect(po))
    @test popt2m == popt
end;

@testset "update_statepar vector unstructured" begin
    u1t = ComponentVector(L = 10.1)
    pt = ComponentVector(k_L = 1.1, k_R = 1/20.1, m = 2.0)
    test_update_statepar_and_get_paropt(ps, u1, p1, popt, u1t, pt)
end
@testset "update_statepar vector structured" begin
    u1t = ComponentVector(a = (a1 = 1, a2 = (a21 = 1, a22 = 2)))
    pt = ComponentVector(b = (b1 = 0.1, b2 = 0.2), c = [1.0, 2.0])
    test_update_statepar_and_get_paropt(psc, u1c, p1c, poptc, u1t, pt)
end
@testset "update_statepar Svector structured" begin
    u1t = ComponentVector(a = (a1 = 1, a2 = (a21 = 1, a22 = 2)))
    pt = ComponentVector(b = (b1 = 0.1, b2 = 0.2), c = [1.0, 2.0])
    test_update_statepar_and_get_paropt(psc, u1s, p1s, poptc, u1t, pt)
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
    p = (p1=1.1,p2=2)
    tspan = (0.0,1.0)
    prob = ODEProblem(f,SVector(Tuple(u0)),tspan,SVector(Tuple(p)))
    #sol = solve(prob)
    #sol[end]
    ps = ProblemParSetterComp1(Axis(keys(u0)),Axis(keys(p)),Axis((:u1,:p2)))
    popt = ComponentVector(u1=1/4, p2=1.2)
    #
    # update_statepar
    prob2 = update_statepar(ps, popt, prob)
    @test prob2.u0[1] == popt.u1
    @test prob2.p[1] == p.p1 # not updated
    @test prob2.p[2] == popt.p2
    #
    # get_paropt
    @test get_paropt(ps, prob2) == getdata(popt) #Tuple(popt) # Vector because provided p as SVector <: AbstractVector
    @test get_paropt_labeled(ps, prob2) == popt 
    @test label_par(ps, prob.p) == ComponentVector(p)
end;

@testset "gradient with Tuple" begin
    fcost = (popt) -> begin
        u0, p = update_statepar(ps, popt, u1, p1)
        d = sum(get_paropt(ps, u0, p))
        d*d
    end
    fcost(popt)
    @test typeof(ForwardDiff.gradient(fcost, popt)) == typeof(popt)
end;

@testset "gradient with SVector" begin
    fcost = (popt) -> begin
        u0, p = update_statepar(ps, popt, u1, p1)
        d = sum(get_paropt(ps, u0, p))
        d*d
    end
    poptsv = SVector(popt)
    fcost(poptsv)
    @test typeof(ForwardDiff.gradient(fcost, poptsv)) == typeof(poptsv)
end;

@testset "gradient with AbstractVector" begin
    fcost = (popt) -> begin
        u0, p = update_statepar(ps, popt, u1, p1)
        d = sum(get_paropt(ps, u0, p))
        d*d
    end
    poptv = collect(popt)
    fcost(poptv)
    @test typeof(ForwardDiff.gradient(fcost, poptv)) == typeof(poptv)
end;

function test_system(ps1, popt_names, m)
    #@infiltrate
    #@code_warntype ProblemParSetterComp1(m, Axis(popt_names))
    #@descend_code_warntype ProblemParSetterComp1(m, Axis(popt_names))
    @test @inferred(symbols_state(ps1)) == (:x, :RHS)
    @test @inferred(symbols_par(ps1)) == (:τ,)
    @test @inferred(symbols_paropt(ps1)) == popt_names
end
@testset "construct from ODESystem" begin
    @named m = samplesystem()
    popt_names = (:RHS, :τ)
    ps1 = ProblemParSetterComp1(m, popt_names)
    # only type stable after function boundary
    test_system(ps1, popt_names, m)
    #
    @named em = embed_system(m, simplify=false)
    ps1e = ProblemParSetterComp1(em, popt_names; strip=true)
    test_system(ps1e, popt_names, em)
end;

@testset "name_paropt" begin
    xn = @inferred name_paropt(ps, collect(1:count_paropt(ps)))
    @test names(xn)[1] == collect(symbols_paropt(ps))
    @test xn[:k_R] == 3 
    xn = @inferred name_state(ps, collect(1:count_state(ps)))
    @test names(xn)[1] == collect(symbols_state(ps))
    xn = @inferred name_par(ps, collect(1:count_par(ps)))
    @test names(xn)[1] == collect(symbols_par(ps))
    #
    frandsym = () -> begin
        syms_arr = rand([:L,:k_L,:k_R],2)
        ntuple(i -> syms_arr[i], 2)
        rand() > 0.5 ? (:L,:k_L,:k_R) : (:L,:k_L)
    end
    #frandsym()
    psr = @inferred ProblemParSetterComp1(keys(u1),keys(p1),frandsym(), Val(false))
    # cannot be inferred, because labels are not known at construction time
    #xl = @inferred label_paropt(psr, collect(1:count_paropt(psr)))
    # but NamedVector is ok
    xn = @inferred name_paropt(psr, collect(1:count_paropt(psr)))
    ftmp = () -> begin
        psr = ProblemParSetterComp1(keys(u1),keys(p1),frandsym(), Val(false))
        xn = name_paropt(psr, collect(1:count_paropt(psr)))
    end
    # @code_warntype ftmp()  # all red!! number of parameters not known
    # need different ProblemParSetterComp1 that does not store integers in type signature
    # and hence will not support LabelledArrays
    _ftmp = (namesopt) -> begin
        local psr = ProblemParSetterComp1(keys(u1),keys(p1),namesopt, Val(false))
        #xn = @inferred label_paropt(psr, collect(1:count_paropt(psr)))
        xn = @inferred name_paropt(psr, collect(1:count_paropt(psr)))
    end
    ftmp = () -> begin
        local namesopt = frandsym()
        xn = @inferred _ftmp(namesopt)
        @show typeof(xn)
        # return value of entire function is not type stable, # because namesopt is not typestable 
        # but inside this function typeof(xn) is known
    end
    xn = ftmp()  

end;








x = ComponentArray(a=1.0, b=[2, 1, 4], c=(a=2, b=[1, 2]))
x[KeepIndex(:c)]

ax = getaxes(x)[1]
ax[:a]

pnames = (:a,:b,:c)
u0names = (:u01, :u02)
poptnames = (:u01, :a, :c)
ax = Axis(poptnames)
CA.indexmap(CA.NullAxis())
popt = [1.0, 2.0, 3.0]
p = [11.0, 12.0, 13.0]
u0 = [21.0, 22.0]
ComponentVector(popt, Axis(poptnames))

pset = ProblemParSetterComp1(u0names, pnames, poptnames)

sym_u0 = (k for (i,k) in enumerate(keys(pset.ax_popt)) if pset.is_u0[i])
sym_p = (k for (i,k) in enumerate(keys(pset.ax_popt)) if !pset.is_u0[i])
[idx(pset.u0k for k in keys(pset.ax_u0)]

tmp = symbols_paropt(pset)
@inferred get_paropt(pset, u0, p)
# get_paropt
g = (pset.is_u0[i] ? u0[idx(pset.ax_u0[k])] : p[idx(pset.ax_p[k])] for (i,k) in enumerate(keys(pset.ax_popt)))
Tuple(g)
Vector(g)
N = count

_ax_symbol = MTKHelpers._ax_symbol
tmp = _ax_symbol(Axis(poptnames))
tmp = _ax_symbol(first(getaxes(ComponentVector(a=(a1=1,a2=2)))))
tmp = _ax_symbol(first(getaxes(ComponentVector(a=(a1=1,a2=(a21=21, a22=22))))))

cv = ComponentVector(a=(a1=1,a2=(a21=21, a22=22)))
tmp = ProblemParSetterComp1(cv, keys(p1), ())
symbols_state(tmp)
axis_state(tmp)


tmp = _ax_symbol(first(getaxes(ComponentVector(a=(a1=1,a2=(a21=21, a22=22))))))
keys(tmp)
t = (a = 1, b = 2)
_extend = (s1,s2) -> isnothing(s2) ? s1 : s1 * "₊" * s2
_extend("a",nothing)
_extend("a","b")
_extend.("a", ["b", nothing, "c"])
_extend.("a", [])


values(popt_names)

cv = ComponentVector(a=1, a=2)
Axis((:a,:a))

ComponentVector([1,2],Axis((:a,:a)))

