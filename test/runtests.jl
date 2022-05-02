using MTKHelpers
using Test

using DifferentialEquations, ModelingToolkit
using StaticArrays, LabelledArrays
using ForwardDiff
#push!(LOAD_PATH, expanduser("~/julia/turingtools/")) # access local package repo
using AxisArrays

#include("test/samplesystem.jl")
include("samplesystem.jl")

@testset "ProblemParSetter" begin
    include("test_problemparsetter.jl")
end

@testset "utilities" begin
    include("test_utils.jl")
end

@testset "smoothstep" begin
    include("test_smoothstep.jl")
end

@testset "cairomakie" begin
    include("test_cairomakie.jl")
end




