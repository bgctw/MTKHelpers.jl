using MTKHelpers
import MTKHelpers as CP
using Test

using DifferentialEquations, ModelingToolkit
using StaticArrays, LabelledArrays
using ForwardDiff
using NamedArrays
# #push!(LOAD_PATH, expanduser("~/julia/turingtools/")) # access local package repo
# import AxisArrays: AxisArray
using ComponentArrays
import ComponentArrays as CA

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


@testset "utilities" begin
    include("test_util.jl")
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




