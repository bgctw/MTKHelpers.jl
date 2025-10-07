using Test
using MTKHelpers
using MTKHelpers: MTKHelpers as CP
using OrdinaryDiffEq, ModelingToolkit
# using ComponentArrays: ComponentArrays as CA
# using StaticArrays: StaticArrays as SA
using NamedArrays: NamedArrays

pkgdir = dirname(dirname(pathof(MTKHelpers)))
include(joinpath(pkgdir,"test","samplesystem.jl"))

@testset "getlast" begin
    @named m = samplesystem()
    @named sys = embed_system(m)
    prob = ODEProblem(sys, vcat([m.x => 0.0], [m.τ => 3.0]), (0.0, 10.0))
    sol = solve(prob, Tsit5())
    res = @inferred getlast(sol, m.x, m.RHS)
    @test res == NamedArrays.NamedArray([sol[m.x, end], sol[m.RHS, end]], ([m.x, m.RHS],))
    res2 = @inferred getlast(sol, res) # providing NamedArray to getlast
    @test res2 == res
end;
