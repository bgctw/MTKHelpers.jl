using Test
using MTKHelpers
using MTKHelpers: MTKHelpers as CP
using OrdinaryDiffEq, ModelingToolkit

using CairoMakie

pkgdir = dirname(dirname(pathof(MTKHelpers)))
include(joinpath(pkgdir,"test","samplesystem.jl"))

@testset "series_sol" begin
    @named m = samplesystem()
    @named me = embed_system(m)
    prob = ODEProblem(me, vcat([m.x => 1.1]), (0.0, 1.0))
    sol = solve(prob, Tsit5())
    fig = Figure()
    ax = CairoMakie.Axis(fig[1, 1]; xlabel = "time (yr)")
    #MTKHelpers.MTKHelpersMakieExt.series_sol!(ax, sol, [m.x, m.RHS])
    series_sol!(ax, sol, [m.x, m.RHS])
    #axislegend(ax, unique=true, position=:lb)
    fig[1, 2] = Legend(fig, ax, unique = true)
    #display(fig)
    @test true
end;

tmp_f = function ()
    ext = Base.get_extension(MTKHelpers, :MTKHelpersMakieExt)
    ext
end
