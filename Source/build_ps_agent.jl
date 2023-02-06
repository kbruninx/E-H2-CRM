function build_ps_agent!(mod::Model)
    # Extract sets
    JH = mod.ext[:sets][:JH]
    JD = mod.ext[:sets][:JD]

    # Extract time series data
    AF = mod.ext[:timeseries][:AF] # avaiability factors

    # Extract parameters
    W = mod.ext[:parameters][:W] # weight of the representative days
    VC = mod.ext[:parameters][:VC] # variable costs
    IC = mod.ext[:parameters][:IC] # annuity investment costs
    if mod.ext[:parameters][:EOM] == 1 && mod.ext[:parameters][:H2] == 1
        η_H2_E = mod.ext[:parameters][:η_H2_E] # efficiency of hydrogen turbines
    end

    # ADMM algorithm parameters
    λ_H2 = mod.ext[:parameters][:λ_H2] # H2 prices
    gH_bar = mod.ext[:parameters][:gH_bar] # element in ADMM penalty term related to hydrogen market
    ρ_H2 = mod.ext[:parameters][:ρ_H2] # rho-value in ADMM related to H2 market
    λ_EOM = mod.ext[:parameters][:λ_EOM] # EOM prices
    g_bar = mod.ext[:parameters][:g_bar] # element in ADMM penalty term related to EOM
    ρ_EOM = mod.ext[:parameters][:ρ_EOM] # rho-value in ADMM related to EOM auctions

    # Create variables
    cap = mod.ext[:variables][:cap] = @variable(mod, lower_bound = 0, base_name = "capacity")
    g = mod.ext[:variables][:g] = @variable(mod, [jh = JH, jd = JD], lower_bound = 0, base_name = "generation")

    # Create affine expressions 
    tot_cost = mod.ext[:expressions][:tot_cost] = @expression(mod,
        IC * cap + sum(W[jd] * VC * g[jh,jd] for jh in JH, jd in JD)
    )

    if mod.ext[:parameters][:EOM] == 1 && mod.ext[:parameters][:H2] == 1
        gH = mod.ext[:expressions][:gH] = @expression(mod,
          - g / η_H2_E
    )
    end

    # Objective 
    if mod.ext[:parameters][:EOM] == 1 && mod.ext[:parameters][:H2] == 1
        mod.ext[:objective] = @objective(mod, Min,
        + tot_cost
        - sum(W[jd] * λ_H2[jh, jd] * gH[jh, jd] for jh in JH, jd in JD)
        - sum(W[jd] * λ_EOM[jh, jd] * g[jh, jd] for jh in JH, jd in JD)
        + sum(ρ_EOM / 2 * W[jd] * (g[jh, jd] - g_bar[jh, jd])^2 for jh in JH, jd in JD)
        + sum(ρ_H2 / 2 * W[jd] * (gH[jh,jd] - gH_bar[jh,jd])^2 for jh in JH, jd in JD) 
    )
    else
        mod.ext[:objective] = @objective(mod, Min,
        + tot_cost
        - sum(W[jd] * λ_EOM[jh, jd] * g[jh, jd] for jh in JH, jd in JD)
        + sum(ρ_EOM / 2 * W[jd] * (g[jh, jd] - g_bar[jh, jd])^2 for jh in JH, jd in JD)
    )
    end   

    # Capacity constraint
    mod.ext[:constraints][:cap_limit] = @constraint(mod, [jh = JH, jd = JD],
        g[jh, jd] <= AF[jh, jd] * cap / 1000  # scaling factor needed to go from GW -> TWh
    )

    return mod
end