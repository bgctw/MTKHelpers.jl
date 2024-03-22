#TestEnv.activate()
using Test
using MTKHelpers
using MTKHelpers: MTKHelpers as CP
using OrdinaryDiffEq, ModelingToolkit
using ComponentArrays: ComponentArrays as CA
using LoggingExtras

using MethodOfLines

#include("testset_utils.jl") # @testset_skip
pkgdir = dirname(dirname(pathof(MTKHelpers)))
include(joinpath(pkgdir, "test", "samplesystem.jl"))

(u0, p, popt, prob, zs, state_pos, pdesys) = get_sys_ex_pde();
pset = ODEProblemParSetter(get_system(prob), popt)

zs
unknowns(get_system(prob)) 
# note that indices 1 and 16 are missing due to boundary conditions
# hence, pset cannot use Symbolics.scalarize to expand Y

@testset "remake" begin
    prob2 = remake(prob, popt, pset)
    @test get_par_labeled(pset, prob2)[keys(popt.par)] == popt.par
    @test get_state_labeled(pset, prob2)[keys(popt.state)] == popt.state
    @test get_paropt_labeled(pset, prob2) == popt
end;
