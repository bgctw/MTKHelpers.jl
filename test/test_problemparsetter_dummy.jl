using MTKHelpers
using MTKHelpers: MTKHelpers as CP

@testset "system with symbolic arrays" begin
    dpset = CP.DummyProblemParSetter()
    @test dpset isa AbstractProblemParSetter
end;

# @testset "error messages on not implemented methods" begin
#     @test_throws ErrorException count_state(dpset)
#     @test_throws ErrorException count_par(dpset)
#     @test_throws ErrorException count_paropt(dpset)
#     #
#     @test_throws ErrorException symbols_state(dpset)
#     @test_throws ErrorException symbols_par(dpset)
#     @test_throws ErrorException symbols_paropt(dpset)
#     #
#     @test_throws ErrorException axis_state(dpset)
#     @test_throws ErrorException axis_par(dpset)
#     @test_throws ErrorException axis_paropt(dpset)
# end;
