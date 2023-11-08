var documenterSearchIndex = {"docs":
[{"location":"system_num_dict/","page":"Translating symbols and Nums","title":"Translating symbols and Nums","text":"CurrentModule = MTKHelpers","category":"page"},{"location":"system_num_dict/#Translating-symbols-and-Nums","page":"Translating symbols and Nums","title":"Translating symbols and Nums","text":"","category":"section"},{"location":"system_num_dict/","page":"Translating symbols and Nums","title":"Translating symbols and Nums","text":"ModelingToolkit constructs and ODEProblem from an ODESystem by supplying Dictionaries that map Nums to values.  However, it is more convenient to store initial states and parameters as ComponentVectors with symbolic keys, instead of Dictionaries with Num keys.","category":"page"},{"location":"system_num_dict/","page":"Translating symbols and Nums","title":"Translating symbols and Nums","text":"The following functions help to translate symbols to Nums of the system.","category":"page"},{"location":"system_num_dict/","page":"Translating symbols and Nums","title":"Translating symbols and Nums","text":"# setting up a simple example composite system and problem\nusing ModelingToolkit, DifferentialEquations, ComponentArrays\nusing MTKHelpers\nfunction samplesystem(;name,τ = 3.0, p1=1.1, p2=1.2) \n    @variables t \n    D = Differential(t) \n    sts = @variables x(t) RHS(t)             # RHS is observed\n    ps = @parameters τ=τ p1=p1 p2=p2       # parameters\n    ODESystem([ RHS  ~ p1/p2 * (1 - x)/τ, D(x) ~ RHS ], t, sts, ps; name)\nend                     \n@named m = samplesystem()\n@named sys = embed_system(m)\n\n# setup initial conditions and parameters as ComponentVectors\nu0 = ComponentVector(m₊x=0.0)\np = ComponentVector(m₊p1=1.11, m₊τ = 3.1,)\n\n# convert to Dict(Num -> value) in order to create the ODEProblem\np_numdict = system_num_dict(p, m)\nprob = ODEProblem(sys, system_num_dict(u0,m), (0,2), p_numdict);\n\npset = ProblemParSetter(sys, ComponentVector()) \np_prob = label_par(pset, prob.p)\np_prob.m₊p2 = 1.2 # from default\np_prob.m₊τ == 3.1 # from p","category":"page"},{"location":"system_num_dict/#system_num_dict","page":"Translating symbols and Nums","title":"system_num_dict","text":"","category":"section"},{"location":"system_num_dict/","page":"Translating symbols and Nums","title":"Translating symbols and Nums","text":"system_num_dict\nget_system_symbol_dict","category":"page"},{"location":"system_num_dict/#MTKHelpers.system_num_dict","page":"Translating symbols and Nums","title":"MTKHelpers.system_num_dict","text":"system_num_dict(d, sys::AbstractSystem, string_sys=nameof(sys))\nsystem_num_dict(d, symbol_dict::AbstractDict)\nsystem_num_dict(d, systems::NTuple)\n\nCreate a Dictionary Num=>value from symbolic Dictionary or ComponentVector.\n\nOmit pairs where no Num was found.\n\nIn the third variant, a tuple of AbstractSystems can be specified to  replace Nums of several (sub-)systems in the dictionary.\n\n\n\n\n\n","category":"function"},{"location":"system_num_dict/#MTKHelpers.get_system_symbol_dict","page":"Translating symbols and Nums","title":"MTKHelpers.get_system_symbol_dict","text":"get_system_symbol_dict(sys::AbstractSystem, string_sys::String=string(nameof(sys)))\nget_system_symbol_dict(systems...)\n\nConstruct a Dict{Symbol => Num} for all properties in sys. All Symbols are prefixed with <string_sys>₊\n\nThe second variant merges the dictionaries obtained from several systems.\n\n\n\n\n\n","category":"function"},{"location":"problemparsetter/","page":"Update parameters","title":"Update parameters","text":"CurrentModule = MTKHelpers","category":"page"},{"location":"problemparsetter/#Setting-Initial-state-and-Parameters-of-a-problem","page":"Update parameters","title":"Setting Initial state and Parameters of a problem","text":"","category":"section"},{"location":"problemparsetter/","page":"Update parameters","title":"Update parameters","text":"Often one wants to change a subset of the initial states,u0, and a subset of parameters,p, of an ODEProblem during an optimization.","category":"page"},{"location":"problemparsetter/","page":"Update parameters","title":"Update parameters","text":"Given u0 and p can be expressed as ComponentVectors,  and popt can be expressed as a ComponentVector of optimized parameters,  which may include initial states, the following class helps updating the corresponding positions in  an ODEProblem.","category":"page"},{"location":"problemparsetter/","page":"Update parameters","title":"Update parameters","text":"Initial states and parameter components must have different names.","category":"page"},{"location":"problemparsetter/","page":"Update parameters","title":"Update parameters","text":"# setting up a simple example composite system and problem\nusing ModelingToolkit, DifferentialEquations, ComponentArrays\nusing MTKHelpers\nfunction samplesystem(;name,τ = 3.0, p1=1.1, p2=1.2) \n    @variables t \n    D = Differential(t) \n    sts = @variables x(t) RHS(t)             # RHS is observed\n    ps = @parameters τ=τ p1=p1 p2=p2       # parameters\n    ODESystem([ RHS  ~ p1/p2 * (1 - x)/τ, D(x) ~ RHS ], t, sts, ps; name)\nend                     \n@named m = samplesystem()\n@named sys = embed_system(m)\nprob = ODEProblem(sys, [m.x => 0.0], (0.0,10.0))\n\n# setup position matching\npopt = ComponentVector(m₊x=0.1, m₊p1=2.1)\npset = ProblemParSetter(sys, popt) # not fully type-inferred\n\n# extract optimized \nget_paropt(pset, prob)          # plain vector\nget_paropt_labeled(pset, prob)  # ComponentVector\nname_paropt(pset, prob)         # NamedVector \n\n# update states and parameters\nprob2 = update_statepar(pset, popt, prob)\nprob2.p # p is still a plain vector\nlabel_par(pset, prob2.p).m₊p1 == popt.m₊p1 # attach labels and access properties\nlabel_state(pset, prob2.u0).m₊x == popt.m₊x # attach labels and access properties\nget_paropt_labeled(pset, prob2) == popt","category":"page"},{"location":"problemparsetter/","page":"Update parameters","title":"Update parameters","text":"Note that ProblemParSetter, pset, is only fully type-inferred when constructed with  three ComponentArrays.Axis objects. This propagates to all ComponentVectors  constructed by it, e.g. with label_state. Hence, its recommended to pass pset across a function barrier for code where performance matters.","category":"page"},{"location":"problemparsetter/#ProblemUpdater","page":"Update parameters","title":"ProblemUpdater","text":"","category":"section"},{"location":"problemparsetter/","page":"Update parameters","title":"Update parameters","text":"A ProblemParSetter can be combined with a KeysProblemParGetter or other specific implementations of AbstractProblemParGetter to  update an ODEProblem based on information already present in the ODEProblem.","category":"page"},{"location":"problemparsetter/","page":"Update parameters","title":"Update parameters","text":"The following example updates parameters k_R and k_P in the ODEProblem to the value of k_L. This can be useful to ensure that these parameters are also changed when optimizing parameter k_L.","category":"page"},{"location":"problemparsetter/","page":"Update parameters","title":"Update parameters","text":"An implementations of AbstractProblemParGetter can use any computation of the source keys to provide its destination keys. It should implement the keys method, so that when constructing the ProblemUpdater, consistent keys are used, as in the example below.","category":"page"},{"location":"problemparsetter/","page":"Update parameters","title":"Update parameters","text":"f = (u,p,t) -> p[1]*u\nu0 = (L=1/2,)\np = (k_L = 1.0, k_R = 2.0, k_P = 3.0)\ntspan = (0.0,1.0)\nprob = ODEProblem(f,collect(u0),tspan,collect(p))\n\nsource = (:k_L,:k_L)\ndest = (:k_R,:k_P)\npu = ProblemUpdater(KeysProblemParGetter(source, dest), keys(u0), keys(p))\nprob2 = pu(prob)\np2 = label_par(par_setter(pu), prob2.p)\np2.k_R == p.k_L\np2.k_P == p.k_L","category":"page"},{"location":"problemparsetter/","page":"Update parameters","title":"Update parameters","text":"ProblemUpdater","category":"page"},{"location":"problemparsetter/#MTKHelpers.ProblemUpdater","page":"Update parameters","title":"MTKHelpers.ProblemUpdater","text":"ProblemUpdater(par_getter, par_setter) \nProblemUpdater(par_getter, u0_keys, p_keys)\n\nEncapsulates updating an ODEProblem based on the problem itself by  Callable (pu::ProblemUpdater)(prob).\n\nMust be initialized with a callable AbstractProblemParGetter,  e.g. KeysProblemParGetter and on a AbstractProblemParSetter, e.g. ProblemParSetter.\n\nThe second form of the constructor creates a ProblemParSetter based on  keys(par_getter).\n\n\n\n\n\n","category":"type"},{"location":"problemparsetter/#ProblemParGetter","page":"Update parameters","title":"ProblemParGetter","text":"","category":"section"},{"location":"problemparsetter/","page":"Update parameters","title":"Update parameters","text":"KeysProblemParGetter","category":"page"},{"location":"problemparsetter/#MTKHelpers.KeysProblemParGetter","page":"Update parameters","title":"MTKHelpers.KeysProblemParGetter","text":"KeysProblemParGetter(source_keys::NTuple{N,Symbol})\n\nProvices callable (pg::KeysProblemParGetter)(pu::ProblemUpdater, prob)].     To be used to get the parameters/state vector to be set by ProblemUpdater.\n\nInitialize with an NTuple of symbols that index into  vcat(label_state(pu.pset, prob.u0), label_par(pu.pset, prob.p)).\n\n\n\n\n\n","category":"type"},{"location":"problemparsetter/","page":"Update parameters","title":"Update parameters","text":"In order to provide computations for parameters to set, declare a concrete subtype of KeysProblemParGetter, and implement a custom  method (pg::MyProblemParGetter)(pu::ProblemUpdater, prob) that returns a vector of parameter values.","category":"page"},{"location":"problemparsetter/#ProblemParSetter","page":"Update parameters","title":"ProblemParSetter","text":"","category":"section"},{"location":"problemparsetter/","page":"Update parameters","title":"Update parameters","text":"The following type stores the necessary information, that can be queried.","category":"page"},{"location":"problemparsetter/","page":"Update parameters","title":"Update parameters","text":"ProblemParSetter\ncount_state\naxis_state\nsymbols_state","category":"page"},{"location":"problemparsetter/#MTKHelpers.ProblemParSetter","page":"Update parameters","title":"MTKHelpers.ProblemParSetter","text":"ProblemParSetter(state_names,par_names,popt_names) \nProblemParSetter(sys::ODESystem, popt_names; strip=false)\n\nHelps keeping track of a subset of initial states and paraemters to be optimized.\n\nArguments\n\nstate_names: ComponentVector or Axis of all the initial states of the problem\npar_names: all the parameters of the problem\npopt_names: the parameter/initial states to be optimized.\n\nIf all of state_names, par_names, and popt_names are type-inferred Axes, then also the constructed ProblemParSetter is type-inferred.\n\nThe states and parameters can be extracted from an ModelingToolkit.ODESystem. If strip=true, then namespaces of parameters of a composed system are removed,  e.g. subcomp₊p becomes p.\n\n\n\n\n\n","category":"type"},{"location":"problemparsetter/#MTKHelpers.count_state","page":"Update parameters","title":"MTKHelpers.count_state","text":"count_state(::AbstractProblemParSetter) \ncount_par(::AbstractProblemParSetter) \ncount_paropt(::AbstractProblemParSetter)\n\nReport the number of problem states, problem parameters and optimized parameters respectively.    \n\n\n\n\n\n","category":"function"},{"location":"problemparsetter/#MTKHelpers.axis_state","page":"Update parameters","title":"MTKHelpers.axis_state","text":"axis_state(pset::AbstractProblemParSetter)\naxis_par(pset::AbstractProblemParSetter)\naxis_paropt(pset::AbstractProblemParSetter)\n\nReport the Axis, i.e. nested component symbols of problem states, problem parameters and  optimized parameters respectively. Returns an AbstractAxis.\n\n\n\n\n\n","category":"function"},{"location":"problemparsetter/#MTKHelpers.symbols_state","page":"Update parameters","title":"MTKHelpers.symbols_state","text":"symbols_state(pset::AbstractProblemParSetter)\nsymbols_par(pset::AbstractProblemParSetter)\nsymbols_paropt(pset::AbstractProblemParSetter)\n\nReport the names, i.e. symbols of problem states, problem parameters and  optimized parameters respectively, i.e. the concatenation of components. Similar to ComponentArrays.label, but inferred from Axis object.   \n\n\n\n\n\n","category":"function"},{"location":"problemparsetter/#Updating-a-problem","page":"Update parameters","title":"Updating a problem","text":"","category":"section"},{"location":"problemparsetter/","page":"Update parameters","title":"Update parameters","text":"update_statepar(::AbstractProblemParSetter,popt, prob::ODEProblem)","category":"page"},{"location":"problemparsetter/#MTKHelpers.update_statepar-Tuple{AbstractProblemParSetter, Any, SciMLBase.ODEProblem}","page":"Update parameters","title":"MTKHelpers.update_statepar","text":"prob_new = update_statepar(ps::ProblemParSetter, popt, prob::ODEProblem) \nu0new, pnew = update_statepar(ps::ProblemParSetter, popt, u0, p)\n\nReturn an updated problem or updates states and parameters where values corresponding to positions in popt hold the values of popt. The type is changed to the promotion type of popt, to allow working with dual numbers.\n\n\n\n\n\n","category":"method"},{"location":"problemparsetter/#Extracting-optimized-parameters","page":"Update parameters","title":"Extracting optimized parameters","text":"","category":"section"},{"location":"problemparsetter/","page":"Update parameters","title":"Update parameters","text":"get_paropt(::AbstractProblemParSetter,  u0, p)","category":"page"},{"location":"problemparsetter/#MTKHelpers.get_paropt-Tuple{AbstractProblemParSetter, Any, Any}","page":"Update parameters","title":"MTKHelpers.get_paropt","text":"get_paropt(pset::AbstractProblemParSetter, prob::ODEProblem; kwargs...)\nget_paropt(pset::AbstractProblemParSetter, u0, p)\n\nget_paropt_labeled(pset::AbstractProblemParSetter, prob::ODEProblem; kwargs...)\nget_paropt_labeled(pset::AbstractProblemParSetter, u0, p)\n\nExtract the initial states and parameters corresponding to the positions that are optimized.     If both u0 and p are AbstractVectors, the result is a Vector, otherwise the result is a Tuple.\n\nThe labeled versions additionally call label_paropt (see label_state)  on the return value.\n\n\n\n\n\n","category":"method"},{"location":"problemparsetter/#Labeling","page":"Update parameters","title":"Labeling","text":"","category":"section"},{"location":"problemparsetter/","page":"Update parameters","title":"Update parameters","text":"label_state(::AbstractProblemParSetter, u0)\nname_state(::AbstractProblemParSetter, state::AbstractVector)","category":"page"},{"location":"problemparsetter/#MTKHelpers.label_state-Tuple{AbstractProblemParSetter, Any}","page":"Update parameters","title":"MTKHelpers.label_state","text":"label_state(pset::AbstractProblemParSetter, u::AbstractVector) \nlabel_par(pset::AbstractProblemParSetter, par::AbstractVector) \nlabel_paropt(pset::AbstractProblemParSetter, popt::AbstractVector)\n\nProduce a labeled version, i.e. a ComponentVector of initial states, parameters, or optimized parameters respectively.\n\n\n\n\n\n","category":"method"},{"location":"problemparsetter/#MTKHelpers.name_state-Tuple{AbstractProblemParSetter, AbstractVector}","page":"Update parameters","title":"MTKHelpers.name_state","text":"name_state(pset, u::AbstractVector) \nname_par(pset, par::AbstractVector) \nname_paropt(pset, popt::AbstractVector) \n\nname_paropt(pset::AbstractProblemParSetter, prob::ODEProblem)\n\nProduce a NamedVector of given state, parameters, or optimized vars.\n\n\n\n\n\n","category":"method"},{"location":"problemparsetter/#Setting-entire-state-and-parameter-vectors","page":"Update parameters","title":"Setting entire state and parameter vectors","text":"","category":"section"},{"location":"problemparsetter/","page":"Update parameters","title":"Update parameters","text":"get_u_map","category":"page"},{"location":"problemparsetter/#MTKHelpers.get_u_map","page":"Update parameters","title":"MTKHelpers.get_u_map","text":"get_u_map(names_u, pset::AbstractProblemParSetter)\nget_p_map(names_p, pset::AbstractProblemParSetter)\n\nMap each state and parameter the ProblemParSetter pset to a position in names.\n\nWhen construction an ODEProblem from a ODESystem, the order of states and  parameters may have changed compared with a previous construction.\n\nIn order to set entire state or parameter vectors, a mapping from current to previous positions, i.e. integer indices, is required,  so that one can get a vectors in the new format by \n\nu0_old[u_map]\np_old[p_map]\n\nThe mapping is constructed by supplying the names of u0old and pold to a ProblemParameterSetter constructed with the current ODESystem.\n\nKeyword arguments\n\ndo_warn_missing: set to true to issue warnings if some ODESystem state or  parameter names are not found in the old names. This may give false warnings for System parameters that have defaults and do not need to be part of the parameter vector.\n\n\n\n\n\n","category":"function"},{"location":"cairomakie/","page":"CairoMakie Helpers","title":"CairoMakie Helpers","text":"CurrentModule = MTKHelpers","category":"page"},{"location":"cairomakie/#CairoMakie-plotting-utils","page":"CairoMakie Helpers","title":"CairoMakie plotting utils","text":"","category":"section"},{"location":"cairomakie/","page":"CairoMakie Helpers","title":"CairoMakie Helpers","text":"series_sol!","category":"page"},{"location":"cairomakie/#MTKHelpers.series_sol!","page":"CairoMakie Helpers","title":"MTKHelpers.series_sol!","text":"series_sol!(ax, sol::AbstractODESolution, vars; tspan=extrema(sol.t), labels=string.(vars), nt=120, kwargs...)\n\ncalls CairoMakie.series for a grid of n points and interpolated values from sol. vars is passed to the solution object and can contain observables. Currently works only with solutions created by a non-composite solver, e.g. Tsit5.\n\n\n\n\n\n","category":"function"},{"location":"embed_system/","page":"Embedding a system","title":"Embedding a system","text":"CurrentModule = MTKHelpers","category":"page"},{"location":"embed_system/#Representing-a-system-as-a-component","page":"Embedding a system","title":"Representing a system as a component","text":"","category":"section"},{"location":"embed_system/","page":"Embedding a system","title":"Embedding a system","text":"embed_system","category":"page"},{"location":"embed_system/#MTKHelpers.embed_system","page":"Embedding a system","title":"MTKHelpers.embed_system","text":"embed_system(m;name)\n\nEmbeds system m as the single component of a larger system. This helps to match the naming of the states, parameters, observables to  the namespace of the system. \n\nusing ModelingToolkit, DifferentialEquations\nusing MTKHelpers\nfunction samplesystem(;name) \n    @variables t \n    D = Differential(t) \n    sts = @variables x(t) RHS(t)  # RHS is observed\n    ps = @parameters τ       # parameters\n    ODESystem([ RHS  ~ (1 - x)/τ, D(x) ~ RHS ], t, sts, ps; name)\nend                     \n@named m = samplesystem()\n@named sys = embed_system(m)\n\n# note that keys are `m.x`,`m.τ` or `m.RHS`.\n# Hence only m needs to be defined rather than all the states, parameters, \n# and observables of the system.\nprob = ODEProblem(sys, [m.x => 0.0], (0.0,10.0), [m.τ => 3.0])\nsol = solve(prob);\ndx1 = sol[m.RHS][1:3] \ndx2 = sol[getproperty(m,:RHS)][1:3]  # access by symbols\n\n# using Plots\n# plot(sol, vars=[m.x,m.RHS])    \n# plot(sol, vars=getproperty.(Ref(m),[:x, :RHS])) # access by symbols\nlength(dx2) == 3\n# output\ntrue\n\n\n\n\n\n","category":"function"},{"location":"embed_system/#Overriding-equations-of-an-existing-system","page":"Embedding a system","title":"Overriding equations of an existing system","text":"","category":"section"},{"location":"embed_system/","page":"Embedding a system","title":"Embedding a system","text":"For debugging bigger systems, it is helful to set some equations to zero or modify/simplify the system in other ways.","category":"page"},{"location":"embed_system/","page":"Embedding a system","title":"Embedding a system","text":"Function override_system takes a set of equations and matches the right-hand site to the equations of the original system and replaces those equations.","category":"page"},{"location":"embed_system/","page":"Embedding a system","title":"Embedding a system","text":"# setting up a simple example composite system and problem\nusing ModelingToolkit, DifferentialEquations, ComponentArrays\nusing MTKHelpers\nfunction samplesystem(;name,τ = 3.0, p1=1.1, p2=1.2) \n    @variables t \n    D = Differential(t) \n    sts = @variables x(t) RHS(t)             # RHS is observed\n    ps = @parameters τ=τ p1=p1 p2=p2       # parameters\n    ODESystem([ RHS  ~ p1/p2 * (1 - x)/τ, D(x) ~ RHS ], t, sts, ps; name)\nend                     \n\n# simplify the system by setting RHS ~ RHS_0 * x\nfunction samplesystem_const(RHS0; name) \n    # demonstrating override_system by setting the RHS to constant first order rate\n    m = samplesystem(;name)\n    @unpack RHS,x = m\n    @parameters t \n    ps = @parameters RHS_0=RHS0\n    D = Differential(t)\n    eqs = [\n        RHS ~ RHS_0 * x,\n    ]\n    sys_ext = override_system(eqs, m; name, ps) \nend  \n\n@named mc = samplesystem_const(-0.1)\n@named sys = embed_system(mc)\nprob = ODEProblem(sys, [mc.x => 1.0], (0.0,10.0))\nsol = solve(prob)\nisapprox(sol(8)[1], exp(-0.1*8), atol = 1e-5)","category":"page"},{"location":"embed_system/","page":"Embedding a system","title":"Embedding a system","text":"override_system","category":"page"},{"location":"embed_system/#MTKHelpers.override_system","page":"Embedding a system","title":"MTKHelpers.override_system","text":"override_system(eqs, basesys::AbstractSystem; \n    name::Symbol=Symbol(string(nameof(basesys))*\"_ext\"), \n    ps=Term[], \n    obs=Equation[], \n    evs=ModelingToolkit.SymbolicContinuousCallback[], \n    defs=Dict()\n)\n\nModify basesys by replacing some equations matched by their left-hand-side. The keyword argument correspond to ODESystem.\n\n\n\n\n\n","category":"function"},{"location":"embed_system/#Utilities","page":"Embedding a system","title":"Utilities","text":"","category":"section"},{"location":"embed_system/","page":"Embedding a system","title":"Embedding a system","text":"symbol\nstrip_namespace\nsymbols_state(::ODESystem)","category":"page"},{"location":"embed_system/#MTKHelpers.symbol","page":"Embedding a system","title":"MTKHelpers.symbol","text":"symbol(t)\n\nExtract the inner symbol from a Term, Num, or BasicSymbolic object.\n\n\n\n\n\n","category":"function"},{"location":"embed_system/#MTKHelpers.strip_namespace","page":"Embedding a system","title":"MTKHelpers.strip_namespace","text":"strip_namespace(s)\n\nOmit the part before the first dot.\n\n\n\n\n\n","category":"function"},{"location":"embed_system/#MTKHelpers.symbols_state-Tuple{ModelingToolkit.ODESystem}","page":"Embedding a system","title":"MTKHelpers.symbols_state","text":"symbols_state(sys::ODESystem)\nsymbols_par(sys::ODESystem)\n\nExtract the basic symbols without namespace of system states and system parameters.\n\n\n\n\n\n","category":"method"},{"location":"zindex/","page":"index","title":"index","text":"CurrentModule = MTKHelpers","category":"page"},{"location":"zindex/#Index","page":"index","title":"Index","text":"","category":"section"},{"location":"zindex/","page":"index","title":"index","text":"","category":"page"},{"location":"solution/","page":"Solution handling","title":"Solution handling","text":"CurrentModule = MTKHelpers","category":"page"},{"location":"solution/#Indexing-into-the-solution","page":"Solution handling","title":"Indexing into the solution","text":"","category":"section"},{"location":"solution/","page":"Solution handling","title":"Solution handling","text":"getlast","category":"page"},{"location":"solution/#MTKHelpers.getlast","page":"Solution handling","title":"MTKHelpers.getlast","text":"getlast(sol, vars...)\ngetlast(sol, vars_vec)\n\nApply getindex for each var to sol and return a NamedArray.\n\nusing ModelingToolkit, DifferentialEquations\nusing MTKHelpers\nusing NamedArrays\nfunction samplesystem(;name) \n    @variables t \n    D = Differential(t) \n    sts = @variables x(t) RHS(t)  # RHS is observed\n    ps = @parameters τ       # parameters\n    ODESystem([ RHS  ~ (1 - x)/τ, D(x) ~ RHS ], t, sts, ps; name)\nend                     \n@named m = samplesystem()\n@named sys = embed_system(m)\nprob = ODEProblem(sys, [m.x => 0.0], (0.0,10.0), [m.τ => 3.0])\nsol = solve(prob);\nres = getlast(sol, m.x, m.RHS)\nres == NamedArray([sol[m.x,end], sol[m.RHS,end]], ([m.x, m.RHS],))   \n# output\ntrue\n\n\n\n\n\n","category":"function"},{"location":"","page":"Home","title":"Home","text":"CurrentModule = MTKHelpers","category":"page"},{"location":"#MTKHelpers","page":"Home","title":"MTKHelpers","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Documentation for MTKHelpers.","category":"page"},{"location":"smoothstep/","page":"Smooth steps","title":"Smooth steps","text":"CurrentModule = MTKHelpers","category":"page"},{"location":"smoothstep/#Smooth-steps","page":"Smooth steps","title":"Smooth steps","text":"","category":"section"},{"location":"smoothstep/","page":"Smooth steps","title":"Smooth steps","text":"ODE solvers have a hard time with step changes, where the derivative changes discontinuously. The following function help to avoid associated problems by approximating the step by a smoother function. Argument dx controls the  smoothness.","category":"page"},{"location":"smoothstep/","page":"Smooth steps","title":"Smooth steps","text":"smoothstep","category":"page"},{"location":"smoothstep/#MTKHelpers.smoothstep","page":"Smooth steps","title":"MTKHelpers.smoothstep","text":"smoothstep(x, x_step, dx, a=zero(x), b=one(x))\n\nsmooth step function: \n\n= a for x <= x_step - dx\n= b for x >= x_step + dx\n= (a-b)/2 for x == x_step\n\nand smooth in between.\n\n\n\n\n\n","category":"function"}]
}
