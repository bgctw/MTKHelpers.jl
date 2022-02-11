using CairoMakie

i_tmp = () -> begin
    @named m = samplesystem()
    @named me = embed_system(m)
    prob = ODEProblem(me, [m.x => 1.1], (0.0,1.0), [])
    sol = solve(prob)
    fig,ax = pdf_figure();
    series_sol!(ax, sol, [m.x, m.RHS])

    ms = structural_simplify(m)
    @parameters t
    @variables x(t), RHS(t)
    prob = ODEProblem(ms, [x => 1.1], (0.0,1.0), [])
    sol = solve(prob)
    fig,ax = pdf_figure();
    series_sol!(ax, sol, [x])

    #using Plots
    Plots.plot(sol, vars=[x, RHS])
    sol(0.3, idxs=[RHS])
end
