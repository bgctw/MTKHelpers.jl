tmpf = () -> begin
    TestEnv.activate()
end

using Test, SafeTestsets
const GROUP = get(ENV, "GROUP", "All") # defined in in CI.yml
@show GROUP

using LoggingExtras
() -> begin
    warn2error_logger = TransformerLogger(global_logger()) do log
        return log.level === Logging.Warn ?  merge(log, (; level=Logging.Error)) : log
    end
    #global_logger(warn2error_logger)
    with_logger(warn2error_logger) do
        #include("test/mytest.jl")
        @safetestset "Tests" include("test/test_odeproblemparsetter.jl")    
    end
end

@time begin
    if GROUP == "All" || GROUP == "Basic"
        #join_path(test_path, ...) does not work, because test_path is unknown in new module
        #@safetestset "Tests" include("test/test_symbolicarray.jl")
        #TODO array parameters @time @safetestset "test_symbolicarray" include("test_symbolicarray.jl")
        #@safetestset "Tests" include("test/test_util_componentarrays.jl")
        #TODO array parameters @time @safetestset "test_util_componentarrays" include("test_util_componentarrays.jl")
        #@safetestset "Tests" include("test/test_problemparsetter_dummy.jl")
        @time @safetestset "test_problemparsetter_dummy" include("test_problemparsetter_dummy.jl")
        #@safetestset "Tests" include("test/test_nullodeproblemparsetter.jl")
        @time @safetestset "test_nullodeproblemparsetter" include("test_nullodeproblemparsetter.jl")
        #@safetestset "Tests" include("test/test_odeproblemparsetter.jl")
        @time @safetestset "test_odeproblemparsetter" include("test_odeproblemparsetter.jl")
        #@safetestset "Tests" include("test/test_odeproblemparsetterconcrete.jl")
        @time @safetestset "test_odeproblemparsetterconcrete" include("test_odeproblemparsetterconcrete.jl")
        #@safetestset "Tests" include("test/test_problemupdater.jl")
        @time @safetestset "test_problemupdater" include("test_problemupdater.jl")
        #@safetestset "Tests" include("test/test_util.jl")
        @time @safetestset "test_util" include("test_util.jl")
        #@safetestset "Tests" include("test/test_prior_util.jl")
        @time @safetestset "test_prior_util" include("test_prior_util.jl")
        #@safetestset "Tests" include("test/test_smoothstep.jl")
        @time @safetestset "test_smoothstep" include("test_smoothstep.jl")
        #@safetestset "Tests" include("test/test_solution.jl")
        @time @safetestset "test_solution" include("test_solution.jl")
        #@safetestset "Tests" include("test/test_util_nums.jl")
        @time @safetestset "test_util_nums" include("test_util_nums.jl")
        #@safetestset "Tests" include("test/test_cvwriter.jl")
        @time @safetestset "test_cvwriter" include("test_cvwriter.jl")
    end

    if GROUP == "All" || GROUP == "PDE"
        #TODO reactivate after fixing bug in MTK/Symbolic/MOL
        # #@safetestset "Tests" include("test/test_util_pde.jl")
        # @time @safetestset "test_util_pde" include("test_util_pde.jl")
        # #@safetestset "Tests" include("test/test_pde.jl")
        # @time @safetestset "test_pde" include("test_pde.jl")
    end

    if GROUP == "All" || GROUP == "Plot"
        #@safetestset "Tests" include("test/test_cairomakie.jl")
        @time @safetestset "test_cairomakie" include("test_cairomakie.jl")
    end

    if GROUP == "All" || GROUP == "JET"
        #@safetestset "Tests" include("test/test_JET.jl")
        #@time @safetestset "test_JET" include("test_JET.jl")
        if VersionNumber("1.11") <= VERSION 
            #@safetestset "test" include("test/test_aqua.jl")
            @time @safetestset "test_aqua" include("test_aqua.jl")
        end
    end
end
