"""
    samplesystem(; name, τ = 3.0, i = 0.1, p1 = 1.1, p2 = 1.2)

Defines a simple ODESystem of exponential decay to p1/p2 with rate τ
"""    
function samplesystem(; name, τ = 3.0, i=0.1, p1 = 1.1, p2 = 1.2)
    @variables t
    D = Differential(t)
    sts = @variables x(t) RHS(t)  # RHS is observed
    ps = @parameters τ=τ p1=p1 p2=p2 i=i       # parameters
    ODESystem([
            RHS ~ i - p1 * x^2 + (p2 - x) / τ,
            D(x) ~ RHS,
        ], t, sts, ps; name)
end

function samplesystem_const(RHS0; name)
    # demonstrating override_system by setting the RHS to constant first order rate
    m = samplesystem(; name)
    @unpack RHS, x = m
    @parameters t
    ps = @parameters RHS_0 = RHS0
    D = Differential(t)
    eqs = [RHS ~ RHS_0 * x]
    sys_ext = override_system(eqs, m; name, ps)
end

function example_pde_system(;name)
    #https://docs.sciml.ai/MethodOfLines/dev/tutorials/brusselator/
    @parameters x y t
    @variables u(..) v(..)
    Dt = Differential(t)
    Dx = Differential(x)
    Dy = Differential(y)
    Dxx = Differential(x)^2
    Dyy = Differential(y)^2
    #
    ∇²(u) = Dxx(u) + Dyy(u)
    #
    brusselator_f(x, y, t) = (((x-0.3)^2 + (y-0.6)^2) <= 0.1^2) * (t >= 1.1) * 5.
    #
    x_min = y_min = t_min = 0.0
    x_max = y_max = 1.0
    t_max = 11.5
    #
    α = 10.
    eq = [Dt(u(x,y,t)) ~ 1. + v(x,y,t)*u(x,y,t)^2 - 4.4*u(x,y,t) + α*∇²(u(x,y,t)) + brusselator_f(x, y, t),
           Dt(v(x,y,t)) ~ 3.4*u(x,y,t) - v(x,y,t)*u(x,y,t)^2 + α*∇²(v(x,y,t))]
    domains = [x ∈ Interval(x_min, x_max),
                  y ∈ Interval(y_min, y_max),
                  t ∈ Interval(t_min, t_max)]
    # Periodic BCs
    u0(x,y) = 22(y*(1-y))^(3/2)
    v0(x,y) = 27(x*(1-x))^(3/2)
    bcs = [u(x,y,0) ~ u0(x,y,0),
           u(0,y,t) ~ u(1,y,t),
           u(x,0,t) ~ u(x,1,t),
    
           v(x,y,0) ~ v0(x,y,0),
           v(0,y,t) ~ v(1,y,t),
           v(x,0,t) ~ v(x,1,t)] 
    #
    @named pdesys = PDESystem(eq,bcs,domains,[x,y,t],[u(x,y,t),v(x,y,t)])
end

