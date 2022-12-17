using Interpolations 
using GLM
using Random, LinearAlgebra
using QuantEcon 
using Parameters   
using ProgressMeter

struct TransitionMatrix
    P::Matrix{Float64}       # 4x4 Aggregate Transition Matrix 
    Pz::Matrix{Float64}      # 2x2 Matrix of Aggregate productivity shock
    Peps_gg::Matrix{Float64} # 2x2 Matrix of Idiosyncratic shock conditional on good state to good state
    Peps_bb::Matrix{Float64} # 2x2 Matrix of Idiosyncratic shock conditional on bad state to bad state
    Peps_gb::Matrix{Float64} # 2x2 Matrix of Idiosyncratic shock conditional on good state to bad state
    Peps_bg::Matrix{Float64} # 2x2 Matrix of Idiosyncratic shock conditional on bad state to good state
end

function Compute_transition_matrix(ug::Real, ub::Real, zg_ave_dur::Real, zb_ave_dur::Real, ug_ave_dur::Real, ub_ave_dur::Real, puu_rel_gb2bb::Real, puu_rel_bg2gg::Real)
    pgg = 1-1/zg_ave_dur # probability of remaining in good state
    pbb = 1-1/zb_ave_dur # probability of remaining in bad state
    pgb = 1-pgg # probability of changing from g to b
    pbg = 1-pbb  # probability of changing from b to g

    p00_gg = 1-1/ug_ave_dur # prob. of 0 to 0 cond. on g to g
    p00_bb = 1-1/ub_ave_dur # prob. of 0 to 0 cond. on b to b
    p01_gg = 1-p00_gg # prob. of 0 to 1 cond. on g to g
    p01_bb = 1-p00_bb # prob. of 0 to 1 cond. on b to b

    p00_gb = puu_rel_gb2bb*p00_bb # prob. of 0 to 0 cond. on g to b
    p00_bg = puu_rel_bg2gg*p00_gg # prob. of 0 to 0 cond. on b to g
    p01_gb = 1 - p00_gb # prob. of 0 to 1 cond. on g to b
    p01_bg = 1 - p00_bg # prob. of 0 to 1 cond. on b to g

    p10_gg = (ug - ug*p00_gg)/(1-ug) # prob. of 1 to 0 cond. on  g to g
    p10_bb = (ub - ub*p00_bb)/(1-ub) # prob. of 1 to 0 cond. on b to b
    p10_gb = (ub - ug*p00_gb)/(1-ug) # prob. of 1 to 0 cond. on g to b
    p10_bg = (ug - ub*p00_bg)/(1-ub) # prob. of 1 to 0 cond on b to g
    p11_gg = 1 - p10_gg # prob. of 1 to 1 cond. on  g to g
    p11_bb = 1 - p10_bb # prob. of 1 to 1 cond. on b to b
    p11_gb = 1 - p10_gb # prob. of 1 to 1 cond. on g to b
    p11_bg = 1 - p10_bg # prob. of 1 to 1 cond on b to g

    P = [pgg*p11_gg pgb*p11_gb pgg*p10_gg pgb*p10_gb;
        pbg*p11_bg pbb*p11_bb pbg*p10_bg pbb*p10_bb;
        pgg*p01_gg pgb*p01_gb pgg*p00_gg pgb*p00_gb;
        pbg*p01_bg pbb*p01_bb pbg*p00_bg pbb*p00_bb]
    
    Pz = [pgg pgb; pbg pbb]
    Peps_gg = [p11_gg p10_gg; p01_gg p00_gg]
    Peps_bb = [p11_bb p10_bb; p01_bb p00_bb]
    Peps_gb = [p11_gb p10_gb; p01_gb p00_gb]
    Peps_bg = [p11_bg p10_bg; p01_bg p00_bg]

    transmat = TransitionMatrix(P, Pz, Peps_gg, Peps_bb, Peps_gb, Peps_bg)
    return transmat
end

function KSParameter(;
    β::AbstractFloat = 0.99, α::AbstractFloat=0.36, δ::Real=0.025, θ::Real=1,
    k_min::Real = 0, k_max::Real = 1000,  k_size::Integer=100, K_min::Real=30,    K_max::Real=50,    K_size::Integer=4,
    z_min::Real=0.99, z_max::Real=1.01,  z_size::Integer=2, eps_min::Real=0.0, eps_max::Real=1.0, eps_size::Integer=2,
    ug::AbstractFloat=0.04, ub::AbstractFloat=0.1, zg_ave_dur::Real=8,       
    zb_ave_dur::Real=8, ug_ave_dur::Real=1.5, ub_ave_dur::Real=2.5, puu_rel_gb2bb::Real=1.25, puu_rel_bg2gg::Real=0.75,
    μ::Real = 0,  update_k::Float64 = 0.7)

    if θ == 1
        utility = LogUtility()
    else
        utility = CRRAUtility(θ)
    end

    l_bar = 1/(1-ub) # Labor
    k_grid = (range(0, stop=k_size-1, length=k_size)/(k_size-1)).^7*(k_max-k_min).+k_min # individual capital grid
    k_grid[1] = k_min; k_grid[end] = k_max; # fit on the k_grid

    K_grid   = range(K_min, stop = K_max, length = K_size)    # aggregate capital grid
    z_grid   = range(z_max, stop = z_min, length = z_size)    # aggregate technology shock
    eps_grid = range(eps_max, stop=eps_min, length=eps_size)  # idiosyncratic employment shock grid
    s_grid   = gridmake(z_grid, eps_grid)                     # shock grid

    transmat = Compute_transition_matrix(ug,ub, zg_ave_dur,zb_ave_dur, ug_ave_dur,ub_ave_dur, puu_rel_gb2bb,puu_rel_bg2gg) # Transition matrices

    ksp = (u = utility, beta = β, alpha = α, delta = δ, theta = θ, l_bar = l_bar, k_min=k_min, k_max=k_max, k_grid=k_grid,
    K_min=K_min, K_max=K_max, K_grid=K_grid, z_grid=z_grid, eps_grid=eps_grid, s_grid=s_grid, k_size=k_size, K_size=K_size, 
    z_size=z_size, eps_size=eps_size, s_size=z_size*eps_size, ug=ug, ub=ub, transmat=transmat, mu=μ, update_k = update_k)
    return ksp
end

# Factor prices, Cobb-Douglas Production Function
r(α::Real, z::Real, K::Real, L::Real)=α*z*K^(α-1)*L^(1-α)
w(α::Real,z::Real,K::Real,L::Real)=(1-α)*z*K^(α)*L^(-α)    

mutable struct KSSolution
    k_opt::Array{Float64,3} # Policy function for individual capital, aggregate capital, shock state
    value::Array{Float64,3} # Value function individual capital, aggregate capital, shock state
    B::Vector{Float64} # Coefficients on approximate aggregate capital law of motion
    R2::Vector{Float64} # Rsquare of approximate aggregate capital law of motion 
end

function KSSolution(ksp::NamedTuple)
    k_opt = ksp.beta*repeat(ksp.k_grid,outer=[1,ksp.K_size,ksp.s_size])
    k_opt = 0.9*repeat(ksp.k_grid,outer=[1,ksp.K_size,ksp.s_size])
    k_opt .= clamp.(k_opt, ksp.k_min, ksp.k_max) # Initial policy function
    
    value = ksp.u.(0.1/0.9*k_opt)/(1-ksp.beta) # Initial value function
    B = [0.0, 1.0, 0.0, 1.0] # Initial value of the coefficients 
    kss = KSSolution(k_opt, value, B, [0.0, 0.0])
    return kss
end

function generate_shocks(ksp::NamedTuple; z_shock_size::Integer = 1100, population::Integer = 10000)

    @unpack Peps_gg, Peps_bg, Peps_gb, Peps_bb = ksp.transmat
    zi_shock = simulate(MarkovChain(ksp.transmat.Pz), z_shock_size) # Simulate Aggregate shocks 
    epsi_shock = Array{Int}(undef, z_shock_size, population) 

    rand_draw = rand(population)     # First period
    if zi_shock[1] == 1 # if good (eps employed)
        epsi_shock[1, :] .= (rand_draw .< ksp.ug) .+ 1 
    elseif zi_shock[1] == 2 # if bad (eps unemployed)
        epsi_shock[1, :] .= (rand_draw .< ksp.ub) .+ 1 
    else
        error("Review z_shocks[1]")
    end

    # From the second period onward   
    for t = 2:z_shock_size 
        draw_eps_shock!(Val(zi_shock[t]), Val(zi_shock[t-1]), view(epsi_shock, t, :), epsi_shock[t-1, :], ksp.transmat)
    end

    # Adjustment of number of employed and unemployed agents 
    for t=1:z_shock_size
        n_e = count(epsi_shock[t,:].==1) # count number of employed
        empl_rate_ideal = ifelse(zi_shock[t] == 1, 1.0-ksp.ug, 1.0-ksp.ub)
        gap = round(Int, empl_rate_ideal*population) - n_e # gap in employment 
        if gap > 0
            become_employed_i = rand(findall(2 .== epsi_shock[t,:]), gap)
            epsi_shock[t, become_employed_i] .= 1
        elseif gap < 0
            become_unemployed_i = rand(findall(1 .== epsi_shock[t, :]), -gap)
            epsi_shock[t, become_unemployed_i] .= 2
        end 
    end
    return zi_shock, epsi_shock    
end

draw_eps_shock!(zi::Val{1}, zi_lag::Val{1}, epsi,  epsi_lag::AbstractVector, transmat::TransitionMatrix) =  draw_eps_shock!(epsi, epsi_lag, transmat.Peps_gg)
draw_eps_shock!(zi::Val{1}, zi_lag::Val{2}, epsi, epsi_lag::AbstractVector, transmat::TransitionMatrix) = draw_eps_shock!(epsi, epsi_lag, transmat.Peps_bg)
draw_eps_shock!(zi::Val{2}, zi_lag::Val{1}, epsi, epsi_lag::AbstractVector, transmat::TransitionMatrix) = draw_eps_shock!(epsi, epsi_lag, transmat.Peps_gb)
draw_eps_shock!(zi::Val{2}, zi_lag::Val{2}, epsi,  epsi_lag::AbstractVector, transmat::TransitionMatrix) = draw_eps_shock!(epsi, epsi_lag, transmat.Peps_bb)

function draw_eps_shock!(epsi_shocks, epsi_shock_before, Peps::AbstractMatrix)
    for i ∈ eachindex(epsi_shocks)
        rand_draw=rand()
        epsi_shocks[i] = ifelse(epsi_shock_before[i] == 1,
                (Peps[1, 1] < rand_draw) + 1,  # if employed before
                (Peps[2, 1] < rand_draw) + 1)  # if unemployed before
    end
    return nothing
end

function solve_ump!(ksp::NamedTuple, kss::KSSolution; max_iter::Integer = 10000, tol::AbstractFloat = 1e-8)
    @unpack alpha, beta, delta, theta, l_bar, mu, update_k, k_grid, k_size, K_grid, K_size,  s_grid, s_size, k_min, k_max = ksp
    global counter = 1
    k_opt_n = similar(kss.k_opt)
    prog = ProgressThresh(tol, "Solving individual UMP by Euler method: ")
    while true
        for s_i = 1:s_size
            z, eps = s_grid[s_i, 1], s_grid[s_i, 2]
            for (K_i, K) = enumerate(K_grid)
                Kp, L = compute_Kp_L(K,s_i,kss.B,ksp)
                for (k_i, k) = enumerate(k_grid)
                    wealth = (r(alpha,z,K,L)+1-delta)*k+w(alpha,z,K,L)*(eps*l_bar + mu*(1-eps)) # compute agents' wealth
                    expec = Expectation_FOC(kss.k_opt[k_i, K_i, s_i], Kp, s_i, ksp)
                    cn = (beta*expec)^(-1.0/theta) # Take the inverse the CRRA utility function 
                    k_opt_n[k_i, K_i, s_i] = wealth - cn # Compute optimal individual capital given cn 
                end
            end
        end
        k_opt_n .= clamp.(k_opt_n, k_min, k_max)
        dif_k = maximum(abs, k_opt_n - kss.k_opt)
        
        ProgressMeter.update!(prog, dif_k)
        if dif_k < tol
            break
        end
        kss.k_opt .= ksp.update_k*k_opt_n .+ (1-ksp.update_k)*kss.k_opt
        global counter += 1
    end
    return counter
end

function Expectation_FOC(kp::Real, Kp::Real, s_i::Integer, ksp::NamedTuple)
    @unpack  alpha, theta, delta, l_bar, mu = ksp 
    @unpack P = ksp.transmat

    global expectation_term = 0.0
    for s_n_i = 1:ksp.s_size
        zp, epsp = ksp.s_grid[s_n_i, 1], ksp.s_grid[s_n_i, 2]
        Kpp, Lp = compute_Kp_L(Kp, s_n_i, kss.B, ksp)
        rn = r(alpha, zp, Kp, Lp)
        kpp = interpolate((ksp.k_grid, ksp.K_grid), kss.k_opt[:, :, s_n_i], Gridded(Linear()))
        cp = (rn+1-delta)*kp + w(alpha, zp ,Kp, Lp)*(epsp*l_bar+mu*(1.0-epsp)) - kpp(kp, Kp)
        global expectation_term +=  P[s_i, s_n_i]*(cp)^(-theta)*(1-delta+rn)
    end 
    return expectation_term
end

function compute_Kp_L(K::Real, s_i::Integer, B::AbstractVector, ksp::NamedTuple)
    Kp, L = ifelse(s_i%ksp.eps_size == 1,
        (exp(B[1]+B[2]*log(K)), ksp.l_bar*(1-ksp.ug)), #  compute aggregate capital and labor 
        (exp(B[3]+B[4]*log(K)), ksp.l_bar*(1-ksp.ub))) 
    Kp = clamp(Kp, ksp.K_min, ksp.K_max)
    return Kp, L
end

struct Stochastic
    epsi_shocks::Matrix{Int}
    k_population::Vector{Float64}
end 
Stochastic(epsi_shocks::Matrix{Int}) =  Stochastic(epsi_shocks, fill(40, size(epsi_shocks, 2))) 

function simulate_aggregate_path!(ksp::NamedTuple, kss::KSSolution, zi_shocks::AbstractVector, K_ts::Vector, sm::Stochastic)
    @unpack epsi_shocks, k_population = sm 
    T = length(zi_shocks)   # simulated duration
    N = size(epsi_shocks, 2) # number of agents

    for (t, z_i) = enumerate(zi_shocks) # Loop over all periods 
        K_ts[t] = mean(k_population) # Aggregate Capital over the distribution of agents
        
        for (i, k) in enumerate(k_population) # Loop over individuals
            eps_i = epsi_shocks[t, i]   # Idiosyncratic shock (agent i, period t)
            s_i = epsi_zi_to_si(eps_i, z_i, ksp.z_size) 
            itp_pol = interpolate((ksp.k_grid, ksp.K_grid), kss.k_opt[:, :, s_i], Gridded(Linear())) # Obtain next period capital
            k_population[i] = itp_pol(k, K_ts[t])
        end
    end
    return nothing
end
epsi_zi_to_si(eps_i::Integer, z_i::Integer, z_size::Integer) = z_i + ksp.z_size*(eps_i-1)

function regression_ALM!(ksp::NamedTuple, kss::KSSolution, zi_shocks::Vector, K_ts::Vector; T_discard::Integer=100) 
    n_g = count(zi_shocks[T_discard+1:end-1] .== 1)
    n_b = count(zi_shocks[T_discard+1:end-1] .== 2)
    
    B_n = Vector{Float64}(undef, 4)

    x_g, y_g = Vector{Float64}(undef, n_g), Vector{Float64}(undef, n_g)
    x_b, y_b = Vector{Float64}(undef, n_b), Vector{Float64}(undef, n_b)    

    global i_g = 1
    global i_b = 1

    # Create the samples for both states 
    for t = T_discard+1:length(zi_shocks)-1
        if zi_shocks[t] == 1
            x_g[i_g] = log(K_ts[t])
            y_g[i_g] = log(K_ts[t+1])
            global i_g += 1
        else
            x_b[i_b] = log(K_ts[t])
            y_b[i_b] = log(K_ts[t+1])
            global i_b += 1
        end
    end
    # Compute regression log(K_ts[t+1]) = β0 + β1 log(K_ts[t]) for both states 
    resg = lm([ones(n_g) x_g], y_g) # Linear regression 
    resb = lm([ones(n_b) x_b], y_b)

    kss.R2 = [r2(resg), r2(resb)]
    B_n[1], B_n[2] = coef(resg) # Store coefficients
    B_n[3], B_n[4] = coef(resb)

    dif_B = maximum(abs, B_n-kss.B)
    return B_n, dif_B
end

function compute_ALM_coef!(sm, ksp::NamedTuple, kss::KSSolution, zi_shocks::Vector{Int}; tol_ump::AbstractFloat=1e-8, max_iter_ump::Integer = 100,
    tol_B::AbstractFloat = 1e-8, max_iter_B::Integer = 20, update_B::AbstractFloat=0.3, T_discard::Integer=100)

    K_ts = Vector{Float64}(undef, length(zi_shocks)) 
    global counter_B = 1
    while true
        println(" --- Iteration over ALM Coefficients: $counter_B ---")
        counter = solve_ump!(ksp, kss, max_iter = max_iter_ump, tol = tol_ump)                    # Solve individual problem
        simulate_aggregate_path!(ksp, kss, zi_shocks, K_ts, sm)                         # Compute aggregate path of capital
        B_n, dif_B = regression_ALM!(ksp, kss, zi_shocks, K_ts, T_discard = T_discard)  # Cbtain new ALM coefficient by regression
 
        if dif_B < tol_B
            println("-----------------------------------------------------")
            println("ALM Coefficients Have Converged : dif = $dif_B")
            println("-----------------------------------------------------")
            break
        end
        kss.B .= update_B .* B_n .+ (1-update_B) .* kss.B # update coefficients 
        global counter_B += 1
    end
    return K_ts, counter_B, counter
end

function plot_Aggregate_Law_Motion_Capital(z_grid::AbstractVector, zi_shocks::Vector, B::Vector, K_ts::Vector; T_discard::Integer = 100)
    
    compute_approxKprime(K, z::Val{1}, B) = exp(B[1]+B[2]*log(K)) # law of motion of capital in good state 
    compute_approxKprime(K, z::Val{2}, B) = exp(B[3]+B[4]*log(K))

    K_ts_approx = similar(K_ts) 
    K_ts_approx[T_discard] = K_ts[T_discard]     # compute approximate ALM for capital

    for t = T_discard:length(zi_shocks)-1
        K_ts_approx[t+1] = compute_approxKprime(K_ts_approx[t], Val(zi_shocks[t]), B)
    end

    p = plot(T_discard+1:length(K_ts), K_ts[T_discard+1:end], lab ="True", color=:red, line=:solid, legend = :topleft)
        plot!(p, T_discard+1:length(K_ts), K_ts_approx[T_discard+1:end],lab="Approximation",color=:blue, line=:dash)
        title!(p, "Aggregate Law of Motion for Capital")
    return p
end

function plot_Evolution_Aggregate_Capital(ksp, kss, K_ts)
    K_min, K_max = minimum(K_ts), maximum(K_ts)
    K_lim = range(K_min, stop=K_max, length=100)
    Kp_g = exp.(kss.B[1] .+ kss.B[2]*log.(K_lim)) # law of motion of capital in good state 
    Kp_b = exp.(kss.B[3] .+ kss.B[4]*log.(K_lim))
    
    p = plot(K_lim, Kp_g, linestyle=:solid, lab="Good State", legend = :topleft)
        plot!(p, K_lim, Kp_b, linestyle=:solid, lab="Bad State")
        plot!(p, K_lim, K_lim, color=:black, linestyle=:dash, lab="45 degree", width=0.5)
        title!(p, "Evolution of Aggregate Capital")
    return p
end
