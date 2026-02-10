using MTKHelpers
using Test
using Aqua

@testset "Code quality (Aqua.jl)" begin
    Aqua.test_all(MTKHelpers;
        unbound_args = false, # does not recognize Union{NTuple{N, Symbol}
        stale_deps = (ignore = [:Requires],),
        ambiguities = false,
        persistent_tasks = false, # fails on CI with PreallocationTools
        )
end;

@testset "ambiguities package" begin
    Aqua.test_ambiguities(MTKHelpers;)
end;
