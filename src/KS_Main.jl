include("KS_Structure.jl")
using Plots     

ksp = KSParameter(); # Create an instance of KSParameter defining parameters of the model
kss = KSSolution(ksp);  # Create an instance of the model solution KSSolution

Random.seed!(0) # set seed 
zi_shocks, epsi_shocks = generate_shocks(ksp;); # Generate aggregate productivity and individual labour shocks

# Solve the model with Euler Method and compute the coefficients of the Aggregate Law of Motion of Capital with linear regression
tol_ump = 1e-1;
tol_B = 1e-1;

sm = Stochastic(epsi_shocks);
# It takes approximately 3 minutes to run the function
K_ts, _, _ = compute_ALM_coef!(sm, ksp, kss, zi_shocks, tol_ump = tol_ump, max_iter_ump = 10000, tol_B = tol_B, max_iter_B = 500, update_B = 0.3, T_discard = 100);

println("Approximated Aggregate Capital Law of Motion by the Agents")
println("The regression log(K_{t+1})= $(round(kss.B[1], digits=4)) + $(round(kss.B[2], digits=4)) log(K_{t}) in good times (R2 = $(round(kss.R2[1], digits=10)))")
println("The regression log(K_{t+1})= $(round(kss.B[3], digits=4)) + $(round(kss.B[4], digits=4)) log(K_{t}) in bad times (R2 = $(round(kss.R2[2], digits=10)))")

# Plots 
plot_Aggregate_Law_Motion_Capital(ksp.z_grid, zi_shocks, kss.B ,K_ts, T_discard = T_discard)
plot_Evolution_Aggregate_Capital(ksp , kss, K_ts)