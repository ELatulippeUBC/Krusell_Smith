module Krusell_Smith

include("KS_Structure.jl")

export KSParameter, KSSolution, generate_shocks, compute_ALM_coef!, 
    plot_Aggregate_Law_Motion_Capital, plot_Evolution_Aggregate_Capital,
    Stochastic, LogUtility, CRRAUtility
end