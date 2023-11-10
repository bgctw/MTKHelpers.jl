function samplesystem(; name, τ = 3.0, p1 = 1.1, p2 = 1.2)
    @variables t
    D = Differential(t)
    sts = @variables x(t) RHS(t)  # RHS is observed
    ps = @parameters τ=τ p1=p1 p2=p2      # parameters
    ODESystem([RHS ~ p1 / p2 * (1 - x) / τ, D(x) ~ RHS], t, sts, ps; name)
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
