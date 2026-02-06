Pkg.activate(;temp=true)
#Pkg.add(["OrdinaryDiffEq","ModelingToolkit","SymbolicIndexingInterface","ComponentArrays","ForwardDiff"])
Pkg.add(["OrdinaryDiffEq","ComponentArrays","ForwardDiff","SciMLStructures"])
#Pkg.develop(url="https://github.com/SciML/SymbolicIndexingInterface.jl")
Pkg.develop("SymbolicIndexingInterface")
Pkg.develop("ModelingToolkit")

using ModelingToolkit
using ModelingToolkit: t_nounits as t, D_nounits as D
using SymbolicIndexingInterface: setp, SymbolicIndexingInterface as SII
using ComponentArrays: ComponentArrays as CA
using ForwardDiff: ForwardDiff
using SciMLStructures
using OrdinaryDiffEq
using SymbolicIndexingInterface


u1 = CA.ComponentVector(L = 10.0)
p1 = CA.ComponentVector(k_L = 1.0, k_R = 1 / 20, m = 2.0)
popt = CA.ComponentVector(par=(m = 4.0,),)
function get_sys1()
    sts = @variables L(t)
    ps = @parameters k_L, k_R, m 
    eq = [D(L) ~ 0, ]
    System(eq, t, sts, vcat(ps...); name=:sys1)
end
sys1 = mtkcompile(get_sys1())

_pmap = Dict( keys(p1) .=> collect(values(p1)) )
prob = ODEProblem(sys1, collect(u1), (0.0,1.1), _pmap)

tmp = (popt) -> begin
  #probo = remake(prob; p = [sys1.m => popt[1]])
  p2 = SciMLStructures.replace(SciMLStructures.Tunable(), prob.p, [sys1.m => popt[1]])
  probo = remake(prob; p = p2)
  d = probo.ps[sys1.m]
  d*d
end
tmp([popt.par.m])
tmp([5.0])
res = ForwardDiff.gradient(tmp, [popt.par.m])


tmpf = () -> begin
  # cannot use setter with ForwardDiff
  setter! = SII.setp(sys1, [sys1.m])
  setter!(prob, [popt.par.m])

  tmp = (popt) -> begin
      setter!(prob, popt)
      d = prob.ps[sys1.m]
      d*d
  end
  tmp([popt.par.m])
  res = ForwardDiff.gradient(tmp, [popt.par.m])
end



prob.ps[:m]
probo.ps[:m]


