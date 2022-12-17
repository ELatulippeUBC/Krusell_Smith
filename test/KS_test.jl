using Test
include("KS_Structure.jl")

println("-----------------------------------------------------")
println("Test on the Parameter Grids")
println("-----------------------------------------------------")
ksp = KSParameter();
kss = KSSolution(ksp); 
@testset "Grid Tests" begin
    @test size(kss.k_opt,1) == length(ksp.k_grid)
    @test size(kss.k_opt,2) == length(ksp.K_grid)
end;

println("-----------------------------------------------------")
println("Test on Functions --- Parameters ")
println("-----------------------------------------------------")

# Test QuantEcon Utility Functions with different paremeters θ
@testset "Grid Tests" begin
    ksp_θ_1 = KSParameter(θ = 1)
    @test ksp_θ_1.u == LogUtility(1.0)

    ksp_θ_2 = KSParameter(θ = 2)
    @test ksp_θ_2.u == CRRAUtility(2.0 )
end;

# Verify that the gridmake function has worked -> It should create a grid of 4 aggregate shocks by 2 individual shocks 
ksp = KSParameter();
@test size(ksp.s_grid) == (4,2)

println("-----------------------------------------------------")
println("Test on Functions --- Maximum Iterations")
println("-----------------------------------------------------")
# Check that the maximum number of iterations was not reached 
max_iter_B = 500 
max_iter_ump = 10000
sm = Stochastic(epsi_shocks);
_, B_counter, counter = compute_ALM_coef!(sm, ksp, kss, zi_shocks, tol_ump = 1e-1, max_iter_ump = 10000, tol_B = 1e-1, max_iter_B = 500, update_B = 0.3, T_discard = 100);
@testset "Maximum Iterations Test" begin
    @test B_counter <= max_iter_B # Regression Method failed to converge with tolerance
    @test counter <= max_iter_ump # Euler Method failed to converge with tolerance
end; 





