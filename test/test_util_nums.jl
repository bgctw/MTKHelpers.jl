@named m = samplesystem()
@named m2 = samplesystem()

@testset "system_num_dict single" begin
    symd = get_system_symbol_dict(m)
    @test symd isa Dict
    @test eltype(keys(symd)) == Symbol
    @test eltype(values(symd)) <:  SymbolicUtils.BasicSymbolic
    @test all((:x, :RHS, :τ) .∈ Ref(keys(symd)))
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
end;

@testset "system_num_dict Tuple" begin
    @parameters t
    sys = compose(ODESystem([], t; name=:sys), [m, m2])
    symd = get_system_symbol_dict(sys)
    @test eltype(keys(symd)) == Symbol
    @test eltype(values(symd)) <:  SymbolicUtils.BasicSymbolic
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
