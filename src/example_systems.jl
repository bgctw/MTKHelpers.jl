function samplesystem_vec(; name, τ = 3.0, i = 0.1, p = [1.1, 1.2, 1.3])
    n_comp = 2
    @parameters t
    D = Differential(t)
    @variables x(..)[1:n_comp] #dx(t)[1:2]  # observed dx now can be accessed
    #sts = @variables x[1:n_comp](t) 
    #ps = @parameters τ=τ p[1:n_comp]=p i=i       # parameters
    ps = @parameters τ=τ i=i p[1:3]=p
    sts = [x(t)[i] for i in 1:n_comp]
    eq = [
        D(x(t)[1]) ~ i - p[1] * x(t)[1] + (p[2] - x(t)[1]^2) / τ,
        D(x(t)[2]) ~ i - p[3] * x(t)[2],
    ]
    #ODESystem(eq, t; name)
    #ODESystem(eq, t, sts, [τ, p[1], p[2], i]; name)
    sys = ODESystem(eq, t, sts, vcat(ps...); name)
    #sys = ODESystem(eq, t, sts, ps; name)
    return sys
end
