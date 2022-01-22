using MTKHelpers
using Test

using DifferentialEquations, ModelingToolkit
using StaticArrays, LabelledArrays
using ForwardDiff

#include("test/samplesystem.jl")
include("samplesystem.jl")

@testset "ProblemParSetter" begin
    include("test_problemparsetter.jl")
end

@testset "utilities" begin
    include("test_utils.jl")
end

