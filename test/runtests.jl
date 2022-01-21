using MTKHelpers
using Test

using DifferentialEquations, ModelingToolkit
using StaticArrays, LabelledArrays
using ForwardDiff

@testset "ProblemParSetter" begin
    include("problemparsetter.jl")
end

@testset "utilities" begin
    include("test_utils.jl")
end

