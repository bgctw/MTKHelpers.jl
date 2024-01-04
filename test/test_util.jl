using Test
using MTKHelpers
using MTKHelpers: MTKHelpers as CP
using OrdinaryDiffEq, ModelingToolkit

test_path = splitpath(pwd())[end] == "test" ? "." : "test"
#include(joinpath(test_path,"samplesystem.jl"))
include("samplesystem.jl")

@testset "symbols_state ODE" begin
    @named m = samplesystem()
    @test symbols_state(m) == [:x, :RHS]
    @test symbols_par(m) == [:τ, :p1, :p2, :i]
end

@testset "simplify_symbol" begin
    @test CP.simplify_symbol("var\"Y_Any[1]\"") == "Y[1]"
    @test CP.simplify_symbol("var\"Y_t_Any[1]\"") == "Y_t[1]"
end;

@testset "symbol_op" begin
    @named m = samplesystem()
    Symbol(m.x)
    @test symbol_op(m.x) == :m₊x # calls Num which calls Par
    @test symbol_op(m.τ) == :m₊τ
    #
    # applying simplify_symbol
    @test symbol_op(Symbol("var\"Y_Any[1]\"")) == Symbol("Y[1]")
end;

@testset "tuplejoin" begin
    @test MTKHelpers.tuplejoin((1, 2), (3,), (4, 5, 6)) == Tuple(i for i in 1:6)
    @test MTKHelpers.tuplejoin((1, 2)) == Tuple(i for i in 1:2)
end;

@testset "strip_namespace" begin
    @test MTKHelpers.strip_namespace(:x) == :x
    @test MTKHelpers.strip_namespace(:m₊x) == :x
    @test MTKHelpers.strip_namespace(:s₊m₊x) == :x
    @test MTKHelpers.strip_namespace("s.m.x") == "x"
end;

# @testset "strip_deriv_num" begin
#     @variables t
#     @variables x(t)
#     Symbol(x)
#     @test strip_deriv_num(x) == :x
# end;

@testset "embed with name different than m" begin
    @named m2 = samplesystem()
    @named sys = embed_system(m2)
    prob = ODEProblem(sys, [m2.x => 0.0], (0.0, 10.0), [m2.τ => 3.0])
    sol = solve(prob, Tsit5())
    @test first(sol[m2.x]) == 0.0
    #plot(sol, vars=[m2.x,m2.RHS])    
    #
    # specify by symbol_op instead of num
    _dict_nums = get_system_symbol_dict(sys)
    prob = ODEProblem(sys, [_dict_nums[:m2₊x] => 0.0], (0.0, 10.0), [m2.τ => 3.0])
end;

@testset "override_system" begin
    @named mc = samplesystem_const(-0.1)
    @named sys = embed_system(mc)
    prob = ODEProblem(sys, [mc.x => 1.0], (0.0, 10.0))
    sol = solve(prob, Tsit5())
    #@test sol(8)[1] ≈ 1-0.1*8
    @test isapprox(sol(8)[1], exp(-0.1 * 8), atol = 1e-5)
end;

@testset "override_system: error on unknown equation" begin
    name = :me
    m = samplesystem(; name)
    @unpack p1, x = m
    @parameters t
    ps = @parameters RHS_0 = -0.1
    D = Differential(t)
    eqs = [
        p1 ~ RHS_0 * x,  # but p1 is not a right-hand-side item
    ]
    @test_throws ErrorException sys_ext=override_system(eqs, m; name, ps)
end;
