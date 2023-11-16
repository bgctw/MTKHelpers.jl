tmp_f = function ()
    pop!(LOAD_PATH)
    push!(LOAD_PATH, expanduser("~/twutz/julia/18_tools/devtools"))
    #
    push!(LOAD_PATH, expanduser("~/twutz/julia/18_tools/scimltools"))
end

using MTKHelpers
import MTKHelpers as CP
using Test
using StaticArrays, LabelledArrays, NamedArrays

using OrdinaryDiffEq, ModelingToolkit
using ForwardDiff
# #push!(LOAD_PATH, expanduser("~/julia/turingtools/")) # access local package repo
# import AxisArrays: AxisArray
using ComponentArrays
using ComponentArrays: ComponentArrays as CA
using Statistics
using Distributions

# https://discourse.julialang.org/t/skipping-a-whole-testset/65006/4
import Test: Test, finish
using Test: DefaultTestSet, Broken
using Test: parse_testset_args
macro testset_skip(args...)
    isempty(args) && error("No arguments to @testset_skip")
    length(args) < 2 && error("First argument to @testset_skip giving reason for "
          *
          "skipping is required")
    skip_reason = args[1]
    desc, testsettype, options = parse_testset_args(args[2:(end - 1)])
    ex = quote
        # record the reason for the skip in the description, and mark the tests as
        # broken, but don't run tests
        local ts = DefaultTestSet(string($desc, " - ", $skip_reason))
        push!(ts.results, Broken(:skipped, "skipped tests"))
        local ret = finish(ts)
        ret
    end
    return ex
end

#include("test/samplesystem.jl")
include("samplesystem.jl")

@testset "symbolicarray" begin
    #include("test/test_symbolicarray.jl")
    include("test_symbolicarray.jl")
end;


@testset "util_componentarrays" begin
    #include("test/test_util_componentarrays.jl")
    include("test_util_componentarrays.jl")
end;

@testset "ProblemParSetter_dummy" begin
    #include("test/test_problemparsetter_dummy.jl")
    include("test_problemparsetter_dummy.jl")
end;

@testset "ODEProblemParSetter" begin
    #include("test/test_problemparsetter.jl")
    include("test_problemparsetter.jl")
end;


@testset "ODEProblemParSetterTyped" begin
    #include("test/test_problemparsettertyped.jl")
    include("test_problemparsettertyped.jl")
end;

@testset "ProblemUpdater" begin
    #include("test/test_problem_updater.jl")
    include("test_problem_updater.jl")
end;

@testset "utilities" begin
    #include("test/test_util.jl")
    include("test_util.jl")
end;

@testset "prior_util" begin
    #include("test/test_prior_util.jl")
    include("test_prior_util.jl")
end;

@testset "smoothstep" begin
    include("test_smoothstep.jl")
end;

@testset "cairomakie" begin
    include("test_cairomakie.jl")
end;

@testset "solution" begin
    include("test_solution.jl")
end;

@testset "util_nums" begin
    include("test_util_nums.jl")
end;

using JET: JET
@testset "JET" begin
    @static if VERSION ≥ v"1.9.2"
        JET.test_package(MTKHelpers; target_modules = (@__MODULE__,))
    end
end;
# JET.report_package(MTKHelpers) # to debug the errors
# JET.report_package(MTKHelpers; target_modules=(@__MODULE__,)) # to debug the errors
