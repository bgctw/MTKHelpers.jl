function samplesystem(;name, τ=3.0, p1=1.1, p2=1.2) 
    @variables t 
    D = Differential(t) 
    sts = @variables x(t) RHS(t)  # RHS is observed
    ps = @parameters τ=τ p1=p1 p2=p2      # parameters
    ODESystem([ RHS  ~ p1/p2 * (1 - x)/τ, D(x) ~ RHS ], t, sts, ps; name)
end       
