tmp_f = function()
    push!(LOAD_PATH, expanduser("~/twutz/julia/makietools"))
    push!(LOAD_PATH, expanduser("~/twutz/julia/18_tools/makietools"))
end
using CairoMakie

@testset "series_sol" begin
    @named m = samplesystem()
    @named me = embed_system(m)
    prob = ODEProblem(me, [m.x => 1.1], (0.0,1.0), [])
    sol = solve(prob,Tsit5())
    fig = Figure(); ax = CairoMakie.Axis(fig[1,1]; xlabel = "time (yr)");
    #MTKHelpers.MTKHelpersMakieExt.series_sol!(ax, sol, [m.x, m.RHS])
    series_sol!(ax, sol, [m.x, m.RHS])
    #axislegend(ax, unique=true, position=:lb)
    fig[1,2] = Legend(fig, ax, unique=true)
    #display(fig)
    @test true
end;

tmp_f = function()
    ext = Base.get_extension(MTKHelpers, :MTKHelpersMakieExt)
    ext
end
