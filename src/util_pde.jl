"""
    grid_exp(n, z_m, efold) 

Generate a Vector of n gridpoints spanning ``[0..z_m]`` with exponentially
increasing distances. The larger the e-foling time, the sparser are the
points near ``z_m``.
"""
function grid_exp(n, z_m, efold=3) 
    grid_exp_i.(1:n, z_m, efold, n)
end
function grid_exp_i(i, z_m, efold, n) 
    x = (i - n) / (1 - n) # transform i in 1..n to x in 1..0
    y = exp(-efold * x)
    # normalize to 0 .. z_m
    (y - exp(-efold)) * 1 / (1 - exp(-efold)) * z_m
end

"""
    Dz_exp(x,x_m,b)
    Iz_exp(x,x_m,b)

Density function decreasing exponentially from zero with e-folding time b 
``x \\in [0,x_m]`` for ``x_m > 0`` for which the intral ``\\int_0^{x_m}  = 1``.
"""
function Dz_exp(x,x_m,b) 
    b/(1-exp(-b*x_m)) * exp(-b*x)
end,
function Iz_exp(x,x_m,b) 
    1/(exp(-b*x_m)-1) * (exp(-b*x) - 1)
end

"""
    Dz_lin(z,z_m)
    Iz_lin(z,z_m)

Constant density for which ``\\int_0^{z_m}  = 1``.
"""
function Dz_lin(z,z_m)
     1/z_m   
end,
function Iz_lin(z,z_m) 
    (z_m-z)/z_m 
end


"""
    get_system(prob::AbstractODEProblem) 

Get the System associated to a problem.
"""
function get_system(prob::AbstractODEProblem)
    :sys ∈ propertynames(prob.f) || error(
        "Cannot get System for a Problem without associated system")
   prob.f.sys
end

"""
    get_discrete_space(prob::AbstractODEProblem) 
    get_discrete_space(sys::AbstractSystem) 

Get the MethodOfLines.DiscreSpace object associated to the system or problem.
"""
function get_discrete_space end

"""
    get_1d_grid(prob::AbstractODEProblem) 
    get_1d_state_grid(sys::AbstractSystem) 

Get spatial grid of a discretized 1-D PDESystem associated with a problem.
The `get_1d_state_grid` variant returns the positions of the state
variables, excluding the boundary conditions.
"""
function get_1d_grid end,
function get_1d_state_grid end

function example_pde_problem end

