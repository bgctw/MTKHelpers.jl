using MethodOfLines
using LoggingExtras

@testset "grid_exp" begin
    grid = grid_exp(5, 10, 4.0) 
    # regression test
    dl = log.(diff(grid))
    #plot(dl)  # straight line
    @test all(isapprox.(diff(dl), first(diff(dl)), rtol=0.01))
end;

prob = LoggingExtras.withlevel(Logging.Error) do
    prob = MTKHelpers.example_pde_problem()
end

@testset "get_system" begin
    sys = get_system(prob);
    @test nameof(sys) == :pdesys2
end;

@testset "get_discrete_space" begin
    ds = get_discrete_space(prob);
    Symbol("Y(t, z)") ∈ Symbol.(ds.ū)
end;

@testset "get_1d_grid" begin
    z_grid = get_1d_grid(prob)
    @test z_grid[1] == 0.0
    @test z_grid[end] == 0.3 #z_m
    z_states = get_1d_state_grid(prob)
    @test z_states == z_grid[2:(end-1)]
end;

@testset "remake problem using ComponentVector" begin
    paropt = ComponentVector(state=[], par=(Y0=201, i_Y_agr=[20.0,30.0]))
    prob2 = remake(prob, paropt)
    #
    axis_p = MTKHelpers.axis_of_nums(Tuple(parameters(get_system(prob))))
    p_old = ComponentVector(prob.p, axis_p)
    p_new = ComponentVector(prob2.p, axis_p)
    @test p_new[keys(paropt.par)] == paropt.par
    otherkeys = Tuple(setdiff(keys(p_old), keys(paropt.par)))
    @test p_new[otherkeys] == p_old[otherkeys]
    @test prob2.u0 == prob.u0
end;

