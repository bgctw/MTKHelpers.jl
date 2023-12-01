@named m = samplesystem()
@named m2 = samplesystem()
@named sys = embed_system(m)

@testset "system_num_dict single" begin
    symd = get_system_symbol_dict(m)
    @test symd isa Dict
    @test eltype(keys(symd)) == Symbol
    @test eltype(values(symd)) <: SymbolicUtils.BasicSymbolic
    @test all((:x, :RHS, :τ) .∈ Ref(keys(symd)))
    cv = ComponentVector(x = 1.2)
    ret = system_num_dict(cv, symd)
    @test all(values(ret) .== [1.2]) # cannot test for x, because only defined in m
    #
    p1 = ComponentVector(x = 1.0, τ = 2.0)
    numd = system_num_dict(p1, m)
    @test numd isa Dict
    num_τ = parameters(m)[1] # defined only inside samplesystem()
    num_x = states(m)[1]
    @test num_τ ∈ keys(numd)
    @test num_x ∈ keys(numd)
    @test numd[num_x] == p1.x
    @test numd[num_τ] == p1.τ
    # 

end;

@testset "system_num_dict Tuple" begin
    @parameters t
    sys = compose(ODESystem([], t; name = :sys), [m, m2])
    symd = get_system_symbol_dict(sys)
    @test eltype(keys(symd)) == Symbol
    @test eltype(values(symd)) <: SymbolicUtils.BasicSymbolic
    @test all((:m₊x, :m₊RHS, :m₊τ) .∈ Ref(keys(symd)))
    @test all((:m2₊x, :m2₊RHS, :m2₊τ) .∈ Ref(keys(symd)))
    #
    p1 = ComponentVector(m₊x = 1.0, m₊τ = 2.0, m2₊x = 3.0, m2₊τ = 4.0)
    numd = system_num_dict(p1, sys)
    @test numd isa Dict
    @test eltype(keys(symd)) == Symbol
    @test all((m.x, m.τ) .∈ Ref(keys(numd)))
    @test all((m2.x, m2.τ) .∈ Ref(keys(numd)))
    @test numd[m.x] == p1.m₊x
    @test numd[m.τ] == p1.m₊τ
    @test numd[m2.x] == p1.m2₊x
    @test numd[m2.τ] == p1.m2₊τ
end;

@testset "base_num" begin
    s = first(states(m))
    ret = base_num(s)
    @test isequal(ret, s)
    sym = :bla
    @test isequal(base_num(sym), sym) # fallback for non-Symbolics
end;

@testset "componentvector_to_numdict" begin
    num_dict_par = get_base_num_dict(parameters(sys))
    cv = ComponentVector(m₊p1 = 10.1)
    ret = MTKHelpers.componentvector_to_numdict(cv, num_dict_par)
    @test ret == Dict(m.p1 => 10.1)
    # test empty subcomponent
    cv2 = ComponentVector(state = [], par = cv)
    ret = MTKHelpers.componentvector_to_numdict(cv2.state, num_dict_par)
    @test ret isa Dict
    @test isempty(ret)
    #prob = ODEProblem(sys, [m.x => 0.0], (0.0,10.0), [m.τ => 3.0])
end;
