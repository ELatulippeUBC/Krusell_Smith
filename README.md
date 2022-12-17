# Krusell Smith Replication - ECON622 UBC 

[![Build Status](https://github.com/ELatulippeUBC/Krusell_Smith_Model.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/ELatulippeUBC/Krusell_Smith_Model.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/ELatulippeUBC/Krusell_Smith_Model.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/ELatulippeUBC/Krusell_Smith_Model.jl)

## Table of contents
* [Information](#Information)
* [Files](#files)
* [Limitations](#limitations)

## Information
This project replicates the model with aggregate uncertainty and idiosyncratic shocks of Krusell-Smith (1998) Income and Wealth Heterogeneity in the Macroeconomy for ECON622. 

The solution strategy is defined as follows:
1. Solve the agent's problem with the Euler Equation Method with interpolation
2. Compute the path of aggregate capital with the policy function found in 1.
3. Update the coefficients of the aggregate law of motion (ALM) of the capital where agents approximate the ALM of capital with a log-linear equation
	$\log(K_{t+1} = \beta_1 + \beta_2 \log (K_{t})$
4. Stop when the coefficients of the true and approximated ALM of capital have converged.

One can improve the approximation of the ALM of capital by reducing the tolerance levels (tol_ump, tol_B) in the src/KS_Main.jl file. 
 
## Files
The project includes the following files: 
* src/KS_Structure.jl : It contains all the main functions to compute the algorithm with a linear regression 
* src/KS_Main.jl : Main file to compute the coefficient values and plots using src/KS_Structure.jl
* src/KS_Structure_NN.jl : It contains all the main functions to compute the algorithm with a neural network
* test/KS_test.jl : It contains all the tests performed on src/KS_Structure.jl

# Limitations
* I could not complete the approximation of the (ALM) of capital with a neural network instead of a log-linear equation. Thus, src/KS_Structure_NN.jl should be improved in future project iterations and only viewed this time as a reference. 
* Even though I followed the procedure in https://julia.quantecon.org/software_engineering/testing.html, I had trouble running the code in the virtual environment since some packages (GLM, Distributions, and others) would not be adequately precompiled. I ran the code on VS Code with Julia 1.6.2
* Moreover, I could not link the repo with Codecov to test my code instead of doing the tests locally only. 
