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

using DifferentialEquations, ModelingToolkit
using ForwardDiff
# #push!(LOAD_PATH, expanduser("~/julia/turingtools/")) # access local package repo
# import AxisArrays: AxisArray
using ComponentArrays
using ComponentArrays: ComponentArrays as CA
using Statistics
using Distributions

#include("test/samplesystem.jl")
include("samplesystem.jl")

@testset "util_componentarrays" begin
    #include("test/test_util_componentarrays.jl")
    include("test_util_componentarrays.jl")
end;

@testset "ProblemParSetter_dummy" begin
    #include("test/test_problemparsetter_dummy.jl")
    include("test_problemparsetter_dummy.jl")
end;

@testset "ProblemParSetter_sym" begin
    #include("test/test_problemparsetter_sym.jl")
    include("test_problemparsetter_sym.jl")
end;

@testset "ProblemParSetter" begin
    #include("test/test_problemparsetter.jl")
    include("test_problemparsetter.jl")
end;

@testset "ProblemUpdater" begin
    #include("test/test_problem_updater.jl")
    include("test_problem_updater.jl")
end;

@testset "utilities" begin
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
