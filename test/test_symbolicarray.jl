
using Test
using MTKHelpers
using MTKHelpers: MTKHelpers as CP
using OrdinaryDiffEq, ModelingToolkit
using ComponentArrays: ComponentArrays as CA
using StaticArrays: StaticArrays as SA

#include("test/testset_utils.jl") # @testset_skip
#include("testset_utils.jl") # @testset_skip

@named m2 = MTKHelpers.samplesystem_vec()
@named sys = embed_system(m2)

@testset "system with symbolic arrays" begin
    st = Symbolics.scalarize(m2.x .=> [1.0, 2.0])
    p_new = Symbolics.scalarize(m2.p .=> [2.1, 2.2, 2.3])
    prob = ODEProblem(sys, st, (0.0, 10.0), p_new)
    @test prob.p == [3.0, 0.1, 2.1, 2.2, 2.3]
    sol = solve(prob, Tsit5())
    # first solution state equals the second entry in the mapping of initial state
    @test first(sol[m2.x]) == last.(st)
    #sol[m2.x[1]]
    #plot(sol, vars=[m2.x,m2.RHS])    
    #
    # specify by symbol_op instead of num
    tmp = @inferred CP.get_scalarized_num_dict(states(sys))
    _dict_nums = get_system_symbol_dict(sys)
    st = Symbolics.scalarize(_dict_nums[:m2₊x] .=> [10.1, 10.2])
    prob = ODEProblem(sys, st, (0.0, 10.0), [_dict_nums[:m2₊τ] => 3.0])
    @test prob.u0 == [10.1, 10.2]
    #
    # specify scalarized x
    pset = ODEProblemParSetter(sys, ())
    u0 = label_state(pset, prob.u0)
    d = system_num_dict(u0, _dict_nums)
    # same keys and contents as in st
    @test length(intersect(keys(d), first.(st))) == length(st)
    @test get.(Ref(d), first.(st), nothing) == last.(st)
    d2 = system_num_dict(CA.ComponentVector(m2₊x= [10.1, 10.2]), _dict_nums)
    @test length(intersect(keys(d2), first.(st))) == length(st)
    @test get.(Ref(d2), first.(st), nothing) == last.(st)
end;

@testset "indices_of_nums and axis_of_nums" begin
    nums = parameters(sys)
    pos_nums = MTKHelpers.indices_of_nums(nums)
    @test pos_nums[3] == (symbol_op(m2.p) => 3:5)
    #
    @test symbol_op(m2.x) == :m2₊x
    @test symbol_op(m2.x[1]) == :m2₊x
    #
    ax = MTKHelpers.axis_of_nums(nums)
    tmp = MTKHelpers.attach_axis((1:MTKHelpers.axis_length(ax)) * 10, ax)
    @test tmp.m2₊p == [30, 40, 50]
    #
    ax = MTKHelpers.axis_of_nums(states(sys))
    tmp = MTKHelpers.attach_axis((1:MTKHelpers.axis_length(ax)) * 10, ax)
    @test tmp.m2₊x == [10, 20]
end;

@testset "validate_keys" begin
    # TODO think about validation and provide system  for creation 
    u1 = CA.ComponentVector(x = 1, y = [1, 2])
    p1 = CA.ComponentVector(a = 1, b = [2, 3], c = 4)
    popt_state = CA.ComponentVector(state = CA.ComponentVector(y = [11, 12]))
    popt_par = CA.ComponentVector(par = CA.ComponentVector(c = 40, b = [12, 13]))
    # valid case, different ordering in par
    pset = get_concrete(ODEProblemParSetter(u1, p1, vcat(popt_state, popt_par)))
    res = @inferred MTKHelpers.validate_keys(pset) #TODO check inferred
    @test res.isvalid
    @test isempty(res.msg)
    # no state keyword 
    @test_throws ErrorException get_concrete(ODEProblemParSetter(u1, p1, popt_par))
    try
        get_concrete(ODEProblemParSetter(u1, p1, popt_par))
    catch e
        @test occursin(r"state", e.msg)
    end
    # no par keyword 
    @test_throws ErrorException get_concrete(ODEProblemParSetter(u1, p1, popt_state))
    try
        get_concrete(ODEProblemParSetter(u1, p1, popt_state))
    catch e
        @test occursin(r"par", e.msg)
    end
    # paropt.state of wrong type
    pset = get_concrete(ODEProblemParSetter(u1, p1,
        vcat(CA.ComponentVector(state = (1:3)), popt_par);
        is_validating = Val(false)))
    res = @inferred MTKHelpers.validate_keys(pset)
    @test !res.isvalid
    @test occursin(r"paropt.state <: ComponentVector", res.msg)
    # paropt.par of wrong type
    pset = get_concrete(ODEProblemParSetter(u1, p1,
        vcat(popt_state, CA.ComponentVector(par = (1:3)));
        is_validating = Val(false)))
    res = @inferred MTKHelpers.validate_keys(pset)
    @test !res.isvalid
    @test occursin(r"paropt.par <: ComponentVector", res.msg)
    # missing state key
    pset = get_concrete(ODEProblemParSetter(u1, p1,
        vcat(CA.ComponentVector(state = CA.ComponentVector(foo = 5)), popt_par);
        is_validating = Val(false)))
    res = @inferred MTKHelpers.validate_keys(pset)
    @test !res.isvalid
    @test occursin(r"part of state", res.msg)
    # wrong length state key
    pset = get_concrete(ODEProblemParSetter(u1, p1,
        vcat(CA.ComponentVector(state = CA.ComponentVector(y = 11)), popt_par);
        is_validating = Val(false)))
    res = @inferred MTKHelpers.validate_keys(pset)
    @test !res.isvalid
    @test occursin(r"length", res.msg)
    # missing par key
    pset = get_concrete(ODEProblemParSetter(u1, p1,
        vcat(popt_state, CA.ComponentVector(par = CA.ComponentVector(foo = 5)));
        is_validating = Val(false)))
    res = @inferred MTKHelpers.validate_keys(pset)
    @test !res.isvalid
    @test occursin(r"part of parameters", res.msg)
    # wrong length par key
    pset = get_concrete(ODEProblemParSetter(u1, p1,
        vcat(popt_state, CA.ComponentVector(par = CA.ComponentVector(c = [41, 42])));
        is_validating = Val(false)))
    res = @inferred MTKHelpers.validate_keys(pset)
    @test !res.isvalid
    @test occursin(r"length", res.msg)
end;

@testset "Parsetter from sys_vec" begin
    #_dict_nums = get_system_symbol_dict(sys)
    st = Symbolics.scalarize(m2.x .=> [1.0, 2.0])
    prob = ODEProblem(sys, st, (0.0, 10.0))
    paropt = CA.ComponentArray(state = (m2₊x = [1.1, 1.2],),
        par = (m2₊p = [10.1, 10.2, 10.3],))
    # test assigning paropt to non-scalarized state and parameters
    paropt_nonsplit = vcat(paropt.state, paropt.par)
    ax_state = MTKHelpers.axis_of_nums(states(sys))
    ax_par = MTKHelpers.axis_of_nums(parameters(sys))
    ax_paropt = first(CA.getaxes(paropt))
    # tmp = attach_axis(collect(1:length(ax_paropt)) * 10, ax_paropt)
    # tmp.state.m2₊x
    #pset1 = get_concrete(ODEProblemParSetter(ax_state, ax_par, paropt))
    pset = ODEProblemParSetter(sys, paropt_nonsplit)
    pset = get_concrete(ODEProblemParSetter(sys, paropt_nonsplit))
    @test axis_paropt(pset) == ax_paropt
    explore_create_SVector = () -> begin
        tmpf = (paropt) -> begin
            # extract as StaticArray
            #SVector{length(paropt)}(paropt) # not inferred
            #SVector{5}(paropt) # inferred
            SVector{MTKHelpers.axis_length(MTKHelpers._get_axis(paropt))}(paropt) # inferred :)
        end
        #@inferred tmpf(paropt)
        CA.ComponentArrays.static_getproperty(paropt.state, Val{:m2₊x}())
    end
    prob_opt = remake(prob, paropt, pset)
    @test typeof(prob_opt.u0) == typeof(prob.u0)
    @test typeof(prob_opt.p) == typeof(prob.p)
    paropt2 = get_paropt_labeled(pset, prob_opt)
    @test paropt2 == paropt
    #
    # initialized pset with symbols
    symbols_state(pset)
    paropt_keys = keys_paropt(pset)
    pset2 = get_concrete(ODEProblemParSetter(sys, paropt_keys))
    @test axis_paropt(pset2) == ax_paropt
end;

@testset "assign_state_par" begin
    u1 = CA.ComponentVector(x = 1, y = [1, 2])
    p1 = CA.ComponentVector(a = 1, b = [2, 3], c = 4)
    paropt = CA.ComponentVector(y = [11, 12], c = 40, b = [12, 13])
    # valid case, different ordering in par
    pset = get_concrete(ODEProblemParSetter(u1, p1, paropt))
    tmp = label_paropt(pset, 1:count_paropt(pset))
    @test tmp.state.y == [1, 2]
    @test tmp.par.c == 3
    @test tmp.par.b == [4, 5]
end

@testset "get_u_map and get_p_map" begin
    u1 = CA.ComponentVector(x = 1.0, y = [2.0, 3.0])
    p1 = CA.ComponentVector(a = 10.0, b = [20.0, 30.0, 40], c = 50)
    pset = get_concrete(ODEProblemParSetter(u1,
        p1,
        CA.ComponentVector(state = u1[CA.KeepIndex(:x)] .* 10,
            par = p1[CA.KeepIndex(:b)] .* 2)))
    # assume that positions have been changed
    u_new = u1[SA.SA[:y, :x]]
    u_map = get_u_map(u_new, pset)
    @test all(u1[u_map] .== u_new)
    #
    # only a subcomponent to be set
    u_new = u1[SA.SA[:x]] .* 2
    u_map = @test_logs (:warn, r":y") get_u_map(u_new, pset; is_warn_missing = true)
    u_up = copy(u1)
    u_up[u_map] .= u_new
    @test all(u_up[u_map] .== u_new)
    missing_keys = setdiff(keys(u1), keys(u_new))
    @test all(u_up[missing_keys] .== u1[missing_keys])
    #
    # specify components by name instead of CA.ComponentVector
    u_new = u1[SA.SA[:y, :x]]
    u_map = get_u_map(keys(u_new), pset)
    @test all(u1[u_map] .== u_new)
    #
    p_new = p1[SA.SA[:b, :a]] .* 2
    p_map = @test_logs (:warn, r":c") get_p_map(p_new, pset; is_warn_missing = true)
    p_up = copy(p1)
    p_up[p_map] .= p_new
    @test all(p_up[p_map] .== p_new)
    missing_keys = setdiff(keys(p1), keys(p_new))
    @test all(p_up[missing_keys] .== p1[missing_keys])
end;
