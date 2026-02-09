#TestEnv.activate()
using Test
using MTKHelpers
using MTKHelpers: MTKHelpers as CP
using OrdinaryDiffEq, ModelingToolkit
using ComponentArrays: ComponentArrays as CA
# using StaticArrays: StaticArrays as SA


pkgdir = dirname(dirname(pathof(MTKHelpers)))
include(joinpath(pkgdir,"test","samplesystem.jl"))
include(joinpath(pkgdir,"test","testset_utils.jl"))

@named m = samplesystem()
@named m2 = samplesystem()
@named sys = embed_system(m)

#@named sys_vec = CP.samplesystem_vec()

@testset "system_num_dict single" begin
    symd = get_system_symbol_dict(m)
    @test symd isa Dict
    @test eltype(keys(symd)) == Symbol
    @test eltype(values(symd)) <: SymbolicUtils.BasicSymbolic
    @test all((:x, :RHS, :τ) .∈ Ref(keys(symd)))
    cv = CA.ComponentVector(x = 1.2)
    ret = @inferred system_num_dict(cv, symd)
    @test all(values(ret) .== [1.2]) # cannot test for x, because only defined in m
    #
    p1 = CA.ComponentVector(x = 1.0, τ = 2.0)
    numd = system_num_dict(p1, m)
    @test numd isa Dict
    num_τ = parameters(m)[1] # defined only inside samplesystem()
    num_x = unknowns(m)[1]
    @test num_τ ∈ keys(numd)
    @test num_x ∈ keys(numd)
    @test numd[num_x] == p1.x
    @test numd[num_τ] == p1.τ
end;

@testset "system_num_dict Tuple" begin
    @parameters t
    sys = compose(System(Equation[], t; name = :sys), [m, m2])
    symd = get_system_symbol_dict(sys)
    @test eltype(keys(symd)) == Symbol
    @test eltype(values(symd)) <: SymbolicUtils.BasicSymbolic
    @test all((:m₊x, :m₊RHS, :m₊τ) .∈ Ref(keys(symd)))
    @test all((:m2₊x, :m2₊RHS, :m2₊τ) .∈ Ref(keys(symd)))
    #
    p1 = CA.ComponentVector(m₊x = 1.0, m₊τ = 2.0, m2₊x = 3.0, m2₊τ = 4.0)
    numd = system_num_dict(p1, sys)
    @test numd isa Dict
    @test eltype(keys(symd)) == Symbol
    @test all((m.x, m.τ) .∈ Ref(keys(numd)))
    @test all((m2.x, m2.τ) .∈ Ref(keys(numd)))
    @test numd[m.x] == p1.m₊x
    @test numd[m.τ] == p1.m₊τ
    @test numd[m2.x] == p1.m2₊x
    @test numd[m2.τ] == p1.m2₊τ
end;

@testset "base_num" begin
    s = first(unknowns(m))
    ret = base_num(s)
    @test isequal(ret, s)
    sym = :bla
    @test isequal(base_num(sym), sym) # fallback for non-Symbolics
end;

@testset "componentvector_to_numdict" begin
    num_dict_par = MTKHelpers.get_base_num_dict(parameters(sys))
    cv = CA.ComponentVector(m₊p1 = 10.1)
    ret = MTKHelpers.componentvector_to_numdict(cv, num_dict_par)
    @test ret == Dict(m.p1 => 10.1)
# end;

# @testset_skip "componentvector_to_numdict empty subcomponent" begin
    # does not work in CA 13.8 (required to solve for Turing)
    cv2 = CA.ComponentVector(state = [], par = cv)
    ret = MTKHelpers.componentvector_to_numdict(cv2.state, num_dict_par)
    @test ret isa Dict
    @test valtype(ret) == Float64
    @test isempty(ret)
    #prob = ODEProblem(sys, [m.x => 0.0], (0.0,10.0), [m.τ => 3.0])
    # 
    ret = MTKHelpers.componentvector_to_numdict(CA.ComponentVector{Float64}(), num_dict_par)
    @test ret isa Dict
    @test valtype(ret) == Float64
    @test isempty(ret)
end;

@testset_skip "expand_base_num_axes" begin
    cv = CA.ComponentVector(p = [2.1, 2.2, 2.3], i = 0.2)
    scalar_num_map = CP.get_scalar_num_map(sys_vec)
    #tmp = scalar_num_map[first(keys(scalar_num_map))]
    cvs = CP.expand_base_num_axes(cv, scalar_num_map)
    @test all(cvs .== cv)
    #@test keys(cvs) == (Symbol("getindex(p, 1)"), Symbol("getindex(p, 2)"),
        #Symbol("getindex(p, 3)"), :i)
    @test keys(cvs) == Symbol.(("p[1]","p[2]","p[3]","i"))
end;

@testset_skip "get_scalarized_num_dict" begin
    nums = unknowns(sys_vec)
    nums2 = vcat(nums, parameters(sys_vec))
    d = CP.get_scalarized_num_dict(nums2)
    @test d[Symbol(nums[1])] === nums[1]
    #
    # see test_symbolicarray for system_num_dict(cv, Dict)
end;
