using MethodOfLines
using LoggingExtras

@testset "grid_exp" begin
    grid = grid_exp(5, 10, 4.0)
    # regression test
    dl = log.(diff(grid))
    #plot(dl)  # straight line
    @test all(isapprox.(diff(dl), first(diff(dl)), rtol = 0.01))
end;

prob = LoggingExtras.withlevel(Logging.Error) do
    prob = MTKHelpers.example_pde_problem()
end

@testset "get_system" begin
    sys = get_system(prob)
    @test nameof(sys) == :pdesys2
end;

@testset "get_discrete_space" begin
    ds = get_discrete_space(prob)
    Symbol("Y(t, z)") ∈ Symbol.(ds.ū)
end;

@testset "get_1d_grid" begin
    z_grid = get_1d_grid(prob)
    @test z_grid[1] == 0.0
    @test z_grid[end] == 0.3 #z_m
    z_states = z_grid[get_1d_state_pos(prob)]
    @test z_states == z_grid[2:(end - 1)]
end;

@testset "remake problem using ComponentVector" begin
    state_pos = get_1d_state_pos(prob)
    Y_new = repeat([200.0], length(state_pos))
    paropt = ComponentVector(state = (Y = Y_new,), par = (Y0 = 201, i_Y_agr = [20.0, 30.0]))
    prob2 = remake(prob, paropt; state_pos)
    #
    axis_p = MTKHelpers.axis_of_nums(Tuple(parameters(get_system(prob))))
    p_old = ComponentVector(prob.p, axis_p)
    p_new = ComponentVector(prob2.p, axis_p)
    @test p_new[keys(paropt.par)] == paropt.par
    otherkeys = Tuple(setdiff(keys(p_old), keys(paropt.par)))
    @test p_new[otherkeys] == p_old[otherkeys]
    @test prob2.u0 == getdata(paropt.state)
end;

@testset "remake problem using ComponentVector u0_array" begin
    proba = LoggingExtras.withlevel(Logging.Error) do
        proba = MTKHelpers.example_pde_problem_arrstate()
    end
    #
    z_grid = get_1d_grid(proba)
    state_pos = get_1d_state_pos(proba)
    n_state = length(state_pos)
    z_m = z_grid[end]
    Y_L = repeat([300.0 / z_m], n_state)
    #
    # need to translate Y[L] to Y[1] key name
    L, R = (1, 2)
    state_names = (Symbol("Y[$i]") for i in (R, L)) # turn around to test position matching
    paropt = ComponentVector(state = (; zip(state_names, [Y_L, Y_L * 3])...),
        par = (Y0 = 201,))
    prob2 = remake(proba, paropt; state_pos)
    #
    axis_p = MTKHelpers.axis_of_nums(Tuple(parameters(get_system(proba))))
    p_old = ComponentVector(proba.p, axis_p)
    p_new = ComponentVector(prob2.p, axis_p)
    @test p_new[keys(paropt.par)] == paropt.par
    otherkeys = Tuple(setdiff(keys(p_old), keys(paropt.par)))
    @test p_new[otherkeys] == p_old[otherkeys]
    axis_u = MTKHelpers.axis_of_nums(Tuple(states(get_system(proba))))
    u_new = ComponentVector(prob2.u0, axis_u)
    @test u_new.var"Y[1]" == getdata(paropt.state.var"Y[1]")
    @test u_new.var"Y[2]" == getdata(paropt.state.var"Y[2]")
end;
