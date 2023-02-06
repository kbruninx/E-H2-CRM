function solve_ps_agent!(m::String,mod::Model)
    # Extract sets
    JH = mod.ext[:sets][:JH]
    JD = mod.ext[:sets][:JD]

    # Extract parameters
    W = mod.ext[:parameters][:W] # weight of the representative days
    VC = mod.ext[:parameters][:VC]
    IC = mod.ext[:parameters][:IC] # overnight investment costs
    λ_EOM = mod.ext[:parameters][:λ_EOM] # EOM prices
    g_bar = mod.ext[:parameters][:g_bar] # element in ADMM penalty term related to EOM
    ρ_EOM = mod.ext[:parameters][:ρ_EOM] # rho-value in ADMM related to EUA auctions
    λ_H2 = mod.ext[:parameters][:λ_H2] # H2 prices
    ρ_H2 = mod.ext[:parameters][:ρ_H2] # rho-value in ADMM related to H2 market
    gH_bar = mod.ext[:parameters][:gH_bar] # element in ADMM penalty term related to hydrogen market

    # Extract variables and expressions
    tot_cost = mod.ext[:expressions][:tot_cost]
    cap = mod.ext[:variables][:cap]
    g = mod.ext[:variables][:g]
    if mod.ext[:parameters][:EOM] == 1 && mod.ext[:parameters][:H2] == 1
        gH = mod.ext[:expressions][:gH]
    end

    # Update objective 
    if m == "electrolysis"
        
        mod.ext[:objective] = @objective(mod, Min,
        + tot_cost
        - sum(W[jd] * λ_H2[jh, jd] * gH[jh, jd] for jh in JH, jd in JD)
        - sum(W[jd] * λ_EOM[jh, jd] * g[jh, jd] for jh in JH, jd in JD)
        + sum(ρ_EOM / 2 * W[jd] * (g[jh, jd] - g_bar[jh, jd])^2 for jh in JH, jd in JD)
        + sum(ρ_H2 / 2 * W[jd] * (gH[jh,jd] - gH_bar[jh,jd])^2 for jh in JH, jd in JD) 
    )
    #elseif for demand agent
    else
        
        mod.ext[:objective] = @objective(mod, Min,
        + tot_cost
        - sum(W[jd] * λ_EOM[jh, jd] * g[jh, jd] for jh in JH, jd in JD)
        + sum(ρ_EOM / 2 * W[jd] * (g[jh, jd] - g_bar[jh, jd])^2 for jh in JH, jd in JD)
    )
    end 

    # solve problem
    optimize!(mod)

    return mod
end