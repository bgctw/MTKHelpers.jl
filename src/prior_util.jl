"""
    fit_Dirichlet_std(p2, σ; σ2 = σ^2, p = vcat(p2, 1 - sum(p2)))

Fit a Dirichlet distribution to (k-1) probabilities p2 so that
standard deviation of a probability 1/k equals σ.

Allow to specify full k-length probability vector p and variance σ2
by keyword arguments to prevent repeated computation.
"""
function fit_Dirichlet_std(p2, σ; σ2 = σ^2, p = vcat(p2, 1 - sum(p2)))
    local n = length(p)
    local s = (n - 1) / (n^2 * σ2) - 1
    Dirichlet(p .* s)
end

"""
    fit_Dirichlet_mode(mode, M)

Fit a Dirichlet distribution so that mode of the distribution correspond to Vector mode.
Need to provide a precision scalar M >= 1. 
Values of M close to 1 result in a flat, while high M result in a peaked distribution.
"""
fit_Dirichlet_mode(mode, M) = mode .* length(mode) .* (M - 1) .+ 1

"""
    simplex_grid(m,n=3)

Create a grid across the n-dimesional unit simplex with m points on each dimension.
If `unit=false` then return an Integer grid ranging from 0 to (m-1).

"""
function simplex_grid(m, n = 3; unit = true)
    n < 2 && error("n needs to be at least 2 but got ", n)
    m < 2 && error("m needs to be at least 2 but got ", m)
    L = num_compositions(m - 1, n)
    A = Array{Int}(undef, L, n)
    simplex_grid!(A, m - 1, n; L)
    unit ? A ./ (m - 1) : A
end

# number of cases for simplex_grid
num_compositions(mstep, n) = binomial(n + mstep - 1, n - 1)

"""
Update matrix A with positions of simplex grid with msteps, i.e. mstep+1 records
for all other dimensions being 0
"""
function simplex_grid!(A, mstep, n; L = num_compositions(mstep, n))
    if n == 2
        A .= cat(0:mstep, mstep:-1:0; dims = 2)
    else
        istart0 = 0
        for im in 0:mstep
            Lsub = num_compositions(mstep - im, n - 1)
            i = istart0 .+ (1:Lsub)
            A[i, 1] .= im
            # update the remaining columns recursively with decreased mstep
            B = @view A[i, 2:n]
            simplex_grid!(B, mstep - im, n - 1; L = Lsub)
            # keep track where to update the next rows
            istart0 = istart0 + Lsub
        end
    end
    nothing
end
