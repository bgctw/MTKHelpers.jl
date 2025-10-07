@testset "fit_Dirichlet_std" begin
    α2 = SA[0.5, 1 / 3]
    σ = 0.2
    d = fit_Dirichlet_std(α2, σ)
    @test mean(d)[1:length(α2)] ≈ α2
    @test sum(mean(d)) ≈ 1.0
    @test var(d)[2] ≈ σ^2
    #using StatsPlots
    #plot(d)
end;

@testset "fit_Dirichlet_mode" begin
    #alphas = fit_Dirichlet_mode(SA[0.6,0.3,0.1], 1.1)
    #alphas = fit_dirichlet_mode(SA[1/3,1/3,1/3], 0.999)
    mode_orig = SA[0.6, 0.3, 0.1]
    alphas = fit_Dirichlet_mode(mode_orig, 1.2)
    d = Dirichlet(alphas)
    @test mode(d) ≈ mode_orig
end;

tmpf = () -> begin
    # plot Dirichlet Distribution
    #using TernaryPlots  # cart2tern
    fd(c_tern; d) = pdf(d, SVector(c_tern))
    function g2(x, y; kwargs...)
        fd(cart2tern(x, y); kwargs...)
    end
    #
    x = range(0, 1, length = 500)
    #grid = reverse(collect(Iterators.product(x,x)), dims=1) # flip for heatmap
    # matrix rows have initial at top, heatmap has 0 at bottom -> ok
    # but change x and y coordinates
    grid = map(x -> SA[x[2], x[1]], collect(Iterators.product(x, x)))
    replace_zero(x) = x == 0 ? missing : x
    # evalue pdf(d) on rectangular grid
    mapg(d) = map(x -> replace_zero(g2(x...; d = d)), grid)
    #A1 = map(x -> x[1], grid); heatmap(A1)
    #A2 = map(x -> x[2], grid); heatmap(A2)
    alphas = fit_Dirichlet_mode(SA[0.6, 0.3, 0.1], 1.2)
    #alphas = fit_dirichlet_mode(SA[1/3,1/3,1/3], 0.999)
    d = Dirichlet(alphas)
    heatmap(mapg(d))
end

@testset "simplex_grid" begin
    A = simplex_grid(2, 2; unit = false)
    @test A == [0 1; 1 0]
    A = simplex_grid(3, 2)
    @test A == transpose([0.0 0.5 1.0; 1.0 0.5 0.0])
    #
    A = simplex_grid(3, 3)
    @test A ==
          [0.0 0.0 1.0; 0.0 0.5 0.5; 0.0 1.0 0.0; 0.5 0.0 0.5; 0.5 0.5 0.0; 1.0 0.0 0.0]
    A = simplex_grid(4, 3; unit = false)
    @test size(A) == (10, 3)
    @test all(sum(A; dims = 2) .== 3)
    #
    A = simplex_grid(4, 4; unit = false)
    @test size(A) == (20, 4)
    @test all(sum(A; dims = 2) .== 3)
end;
