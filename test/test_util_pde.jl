using Test
using MTKHelpers
using MTKHelpers: MTKHelpers as CP
using OrdinaryDiffEq, ModelingToolkit
using ComponentArrays: ComponentArrays as CA
using LoggingExtras

using MethodOfLines

#include("testset_utils.jl") # @testset_skip

@testset "grid_exp" begin
    grid = grid_exp(5, 10, 4.0)
    # regression test
    dl = log.(diff(grid))
    #plot(dl)  # straight line
    @test all(isapprox.(diff(dl), first(diff(dl)), rtol = 0.01))
end;

@testset "Dz_exp" begin
    @test Dz_exp(0.3, 0.3, 2)≈2.43 rtol=0.01 # regression test
    @test Iz_exp(0.3, 0.3, 2) == 1.0
    #
    @test Dz_lin(0.3, 0.3) ≈ 1 / 0.3
    @test Iz_lin(0.3, 0.3) == 1.0
end;

prob = LoggingExtras.withlevel(Logging.Error) do
    prob = MTKHelpers.example_pde_problem()
end

tmpf = () -> begin
    prob.tspan
    sol = solve(prob, Tsit5())
    X = getproperty(get_system(prob), :X; namespace = false)
    sol[X]
end

@testset "get_system" begin
    sys = get_system(prob)
    @test nameof(sys) == :pdesys2
end;

@testset "get_discrete_space" begin
    ds = get_discrete_space(get_system(prob))
    Symbol("Y(t, z)") ∈ Symbol.(ds.ū)
end;

@testset "get_1d_grid" begin
    sys = get_system(prob)
    z_grid = get_1d_grid(sys)
    @test z_grid[1] == 0.0
    @test z_grid[end] == 0.3 #z_m
    z_states = z_grid[get_1d_state_pos(sys)]
    @test z_states == z_grid[2:(end - 1)]
end;

tmpf = () -> begin
    k = last(keys(cv))
    simplify_symbol(k)
    symbol_op(k)
    (; zip(simplify_symbol.(keys(cv)), values(cv))...)
end

@testset "remake problem using CA.ComponentVector" begin
    sys = get_system(prob)
    state_pos = get_1d_state_pos(sys)
    Y_new = 200.0 .+ (1:length(state_pos)) ./ 100
    paropt = CA.ComponentVector(state = (Y = Y_new,),
        par = (Y0 = 201, i_Y_agr = [20.0, 30.0]))
    #
    pset = ODEProblemParSetter(sys, paropt)
    label_state(pset, prob.u0)
    label_par(pset, prob.p)
    label_paropt(pset, CA.getdata(paropt))
    @test :X ∈ symbols_state(pset)   # scalar variables simplified
    #
    prob2 = remake(prob, paropt, pset)
    @test label_par(pset, prob2.p).Y0 == paropt.par.Y0
    # test the order
    @test label_state(pset, prob2.u0)[CA.KeepIndex(1)][1] == paropt.state[:Y][1]
    @test label_state(pset, prob2.u0)[CA.KeepIndex(length(Y_new))][1] ==
          paropt.state[:Y][end]
    #
    paropt2 = get_paropt_labeled(pset, prob2)
    @test paropt2 == paropt
    #prob2 = MTKHelpers.remake_cv(prob, paropt; state_pos)
    #
    # axis_p = MTKHelpers.axis_of_nums(Tuple(parameters(get_system(prob))))
    # p_old = CA.ComponentVector(prob.p, axis_p)
    # p_new = CA.ComponentVector(prob2.p, axis_p)
    # @test p_new[keys(paropt.par)] == paropt.par
    # otherkeys = Tuple(setdiff(keys(p_old), keys(paropt.par)))
    # @test p_new[otherkeys] == p_old[otherkeys]
    # @test prob2.u0 == CA.getdata(paropt.state)
end;

@testset "remake problem using CA.ComponentVector u0_array" begin
    proba = LoggingExtras.withlevel(Logging.Error) do
        proba = MTKHelpers.example_pde_problem_arrstate()
    end
    #
    sys = get_system(proba)
    z_grid = get_1d_grid(sys)
    state_pos = get_1d_state_pos(sys)
    n_state = length(state_pos)
    z_m = z_grid[end]
    Y_R0 = 300.0 / z_m
    Y_R = repeat([Y_R0], n_state) .+ sort(state_pos) ./ 10
    Y_L = Y_R * 3
    #
    # need to translate Y[L] to Y[1] key name
    L, R = (1, 2)
    state_names = (Symbol("Y[$i]") for i in (R, L)) # turn around to test position matching
    paropt = CA.ComponentVector(state = (; zip(state_names, [Y_L, Y_R])...),
        par = (Y0 = 201,))
    #
    pset = ODEProblemParSetter(sys, paropt)
    label_state(pset, proba.u0)
    label_par(pset, proba.p)
    label_paropt(pset, CA.getdata(paropt))
    #
    prob2 = remake(proba, paropt, pset)
    @test label_par(pset, prob2.p).Y0 == paropt.par.Y0
    label_state(pset, prob2.u0)
    # relies on first state pertains to Y[1] == Y[R] - this can change
    #@test label_state(pset, prob2.u0)[CA.KeepIndex(1)][1] == Y_R0 + state_pos[1]/10
    #
    # test only reverse translating positions
    paropt2 = get_paropt_labeled(pset, prob2)
    @test paropt2 == paropt
    #
    # prob2 = MTKHelpers.remake_cv(proba, paropt; state_pos)
    #
    # axis_p = MTKHelpers.axis_of_nums(Tuple(parameters(get_system(proba))))
    # p_old = CA.ComponentVector(proba.p, axis_p)
    # p_new = CA.ComponentVector(prob2.p, axis_p)
    # @test p_new[keys(paropt.par)] == paropt.par
    # otherkeys = Tuple(setdiff(keys(p_old), keys(paropt.par)))
    # @test p_new[otherkeys] == p_old[otherkeys]
    # #
    # tmp = unknowns(get_system(proba))
    # axis_u = MTKHelpers.axis_of_nums(Tuple(unknowns(get_system(proba))))
    # u_new = CA.ComponentVector(prob2.u0, axis_u)
    # @test u_new.var"Y[1]" == CA.getdata(paropt.state.var"Y[1]")
    # @test u_new.var"Y[2]" == CA.getdata(paropt.state.var"Y[2]")
end;

tmp_f = () -> begin
    # setting state by Dict num -> value
    sys = get_system(proba)
    sd = get_system_symbol_dict(sys)
    tmp2 = Symbolics.scalarize(sd[Symbol("Y[1]")])[state_pos] .=>
        identity.(z_grid[state_pos])
    prob3 = remake(proba, u0 = tmp2)
    #
    # translating ComponentVector(vector) to CompoenentVector(scalars...)
    tmp = CP.expand_base_num(sd[Symbol("Y[1]")], state_pos)
    tmp = CP.expand_base_num(sd[Symbol("Y[1]")], sys)
    ax_scalar = CA.Axis(Symbol.(tmp)...)
    sys = get_system(proba)
    s1 = unknowns(sys)[1]
    nums = Tuple(unknowns(get_system(proba)))
    CP._get_axis(Symbol.(nums))
    u_popt = CA.ComponentVector(var"Y[1]" = z_grid[state_pos])
    u_popt_s = CP.attach_axis(u_popt, ax_scalar)
    u_popt_s[Symbol(tmp[1])] == z_grid[state_pos[1]]
    prob3 = remake(proba, u0 = u_popt_s)
    keys(u_popt_s)[1]
    pos_nums = CP.indices_of_nums()
    #
    cv = CA.ComponentVector(var"Y[1]" = z_grid[state_pos], k_Y = 2.1)
    cvs = CP.expand_base_num_axes(cv, sys)
    @test all(cvs .== cv)
    #
    st = unknowns(sys)
    t = NamedTuple(t...)
    cvs = CA.ComponentVector(; zip(Symbol.(st), 1:length(st))...)
    cvs[Symbol.(tmp)]
    keys(cvs)
    # 
    paropt = CA.ComponentVector(state = (var"Y[1]" = z_grid[state_pos],),
        par = (k_Y = 2.1,))
end
