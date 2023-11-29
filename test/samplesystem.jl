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

# function samplesystem_vec(; name, τ = 3.0, i=0.1, p = [1.1, 1.2])
#     ncomp = 2
#     @variables t
#     D = Differential(t)
#     sts = @variables x(t)[1:2] dx(t)[1:2]  # observed dx now can be accessed
#     ps = @parameters τ=τ p p2=p[2] i=i       # parameters
#     ODESystem([
#             dx₁(t) ~ i - p1 * x₁(t)^2 + (p2 - x₁(t)) / τ,
#             dx₂(t) ~ i - p1 * x₂(t)^2 + (p2 - x₂(t)) / τ,
#             D(x₁) ~ dx₁(t),
#             D(x₂) ~ dx₂(t),
#         ], t, sts, ps; name)
# end




