@testset "construct from ODESystem" begin
  @named m = samplesystem()
  #
  symd = get_system_symbol_dict(m)
  @test symd isa Dict{Symbol, Num}
  @test all((:m₊x,:m₊RHS,:m₊τ) .∈ Ref(keys(symd)))
  #
  p1 = ComponentVector(m₊x = 1.0, m₊τ = 2.0)
  numd = system_num_dict(p1, m)
  @test numd isa Dict{Num}
  @test all((m.x,m.τ) .∈ Ref(keys(numd)))
  @test numd[m.x] == p1.m₊x
  @test numd[m.τ] == p1.m₊τ
end;

