function solve_ps_agent!(mod::Model)
    # Extract sets
    JH = mod.ext[:sets][:JH]
    JD = mod.ext[:sets][:JD]
   
    # Extract parameters
    W = mod.ext[:parameters][:W] # weight of the representative days
    VC  = mod.ext[:parameters][:VC]  
    IC = mod.ext[:parameters][:IC] # overnight investment costs
    λ_EOM = mod.ext[:parameters][:λ_EOM] # EOM prices
    g_bar = mod.ext[:parameters][:g_bar] # element in ADMM penalty term related to EOM
    ρ_EOM = mod.ext[:parameters][:ρ_EOM] # rho-value in ADMM related to EUA auctions
  
    # Extract variables and expressions
    cap = mod.ext[:variables][:cap]  
    g = mod.ext[:variables][:g]  

    # Update objective 
    mod.ext[:objective] = @objective(mod, Min,
        + IC*cap
        + sum(W[jd]*VC*g[jh,jd] for jh in JH, jd in JD)
        - sum(W[jd]*λ_EOM[jh,jd]*g[jh,jd] for jh in JH, jd in JD)
        + sum(ρ_EOM/2*W[jd]*(g[jh,jd] - g_bar[jh,jd])^2 for jh in JH, jd in JD)
    )

    # solve problem
    optimize!(mod);

    return mod
end