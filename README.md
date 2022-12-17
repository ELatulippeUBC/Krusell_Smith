# Krusell_Smith_Model

[![Build Status](https://github.com/ELatulippeUBC/Krusell_Smith_Model.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/ELatulippeUBC/Krusell_Smith_Model.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/ELatulippeUBC/Krusell_Smith_Model.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/ELatulippeUBC/Krusell_Smith_Model.jl)

## Table of contents
* [General info](#general-info)
* [Files](#files)

## General info
This project replicates the model with aggregate uncertainty and idiosyncratic shocks of Krusell-Smith (1998) Income and Wealth Heterogeneity in the Macroeconomy for ECON622. 

Solution Strategy is defined as follows:
1. Solve the agent's problem with the Euler Equation Method with interpolation
2. Compute the path of aggregate capital with the policy function found in 1.
3. Update the coefficients of the aggregate law of motion (ALM) of capital where agents approximate the ALM with a log-linear equation
	$\log(K_{t+1} = \beta_1 + \beta_2 \log (K_{t})$
4. Stop when the convergence of the true law of motion and the one approximated by the agents has been reached. 
	
## Files
The project includes the following files: 
* src/KS_Structure.jl : It contains all the main functions to compute the algorithm with a linear regression 
* src/KS_Main.jl : Main file to compute the coefficient values and plots using src/KS_Structure.jl
* src/KS_Structure_NN.jl : It contains all the main functions to compute the algorithm with a neural network
* test/KS_test.jl : It contains all the tests performed on src/KS_Structure.jl

# Limitations
* I was not able to make the neural network functions work to 
* 
	
