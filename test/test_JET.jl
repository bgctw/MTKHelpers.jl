using JET: JET
using MTKHelpers
# using MethodOfLines
# using CairoMakie

@testset "JET" begin
    @static if VERSION ≥ v"1.9.2"
        JET.test_package(MTKHelpers; target_modules = (@__MODULE__,))
    end
end;

tmpf = () -> begin
    JET.report_package(MTKHelpers) # to debug the errors
    JET.report_package(MTKHelpers; target_modules = (@__MODULE__,)) # to debug the errors
end
