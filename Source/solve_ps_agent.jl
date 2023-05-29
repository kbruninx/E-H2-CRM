function solve_ps_agent!(m::String, mod::Model, EOM::Dict)
    # Extract sets
    JH = mod.ext[:sets][:JH]
    JD = mod.ext[:sets][:JD]
    JY = mod.ext[:sets][:JY]

    # Extract parameters
    W = mod.ext[:parameters][:W] # weight of the representative days
    P = mod.ext[:parameters][:P] # probability of each scenario
    γ = mod.ext[:parameters][:γ] # weight of expected revenues and CVAR
    β = mod.ext[:parameters][:β] # risk aversion parametrization
    σ = mod.ext[:parameters][:σ] # switch capacity market
    σH = mod.ext[:parameters][:σH] # switch H2 capacity market

    # ADMM algorithm parameters
    λ_EOM = mod.ext[:parameters][:λ_EOM] # EOM prices
    g_bar = mod.ext[:parameters][:g_bar] # element in ADMM penalty term related to EOM
    ρ_EOM = mod.ext[:parameters][:ρ_EOM] # rho-value in ADMM related to EUA auctions
    λ_H2 = mod.ext[:parameters][:λ_H2] # H2 prices
    ρ_H2 = mod.ext[:parameters][:ρ_H2] # rho-value in ADMM related to H2 market
    gH_bar = mod.ext[:parameters][:gH_bar] # element in ADMM penalty term related to hydrogen market
    λ_CM = mod.ext[:parameters][:λ_CM] # capacity market prices
    cap_bar = mod.ext[:parameters][:cap_bar]  # element in ADMM penalty term related to CM
    ρ_CM = mod.ext[:parameters][:ρ_CM]  # rho-value in ADMM related to CM auctions
    λ_HCM = mod.ext[:parameters][:λ_HCM] # hydrogen capacity market prices
    capH_bar = mod.ext[:parameters][:capH_bar]  # element in ADMM penalty term related to hydrogen CM
    ρ_HCM = mod.ext[:parameters][:ρ_HCM]  # rho-value in ADMM related to hydrogen CM auctions


    # Extract variables and update objective

    if m == "Edemand"
        # Extract parameters
        WTP = mod.ext[:parameters][:WTP]  # willingness to pay ("price cap") of consumers
        ela = EOM["elasticity"]

        # Decision variables
        g_VOLL = mod.ext[:variables][:g_VOLL]
        g_ela = mod.ext[:variables][:g_ela]
        α = mod.ext[:variables][:α]
        u = mod.ext[:variables][:u]

        # Expressions 

        g_positive = mod.ext[:expressions][:g_positive] = @expression(mod, g_VOLL + g_ela)
        g =  mod.ext[:expressions][:g] = @expression(mod, - g_positive)

        profit = mod.ext[:expressions][:profit] = @expression(mod, [jy = JY],      # the welfare is formulated negative to make the problem convex
        sum(W[jd, jy] * (WTP * g_positive[jh, jd, jy] - (g_ela[jh, jd, jy])^2 * WTP / (2*ela[jh, jd, jy]) 
        - λ_EOM[jh, jd, jy] * g_positive[jh, jd, jy]) for jh in JH, jd in JD)
        - σ * λ_CM * CM["D"] )

        cost = mod.ext[:expressions][:cost] = @expression(mod, [jy = JY],
        sum(W[jd, jy] * λ_EOM[jh, jd, jy] * g_positive[jh, jd, jy] for jh in JH, jd in JD))

        CVAR = mod.ext[:expressions][:CVAR] = @expression(mod,
        α - ((1/β) * sum(P[jy] * u[jy] for jy in JY)) )

        # Update objective function 
        mod.ext[:objective] = @objective(mod, Min,
        - γ * sum(P[jy] * profit[jy] for jy in JY) - (1 - γ) * CVAR
        + sum(ρ_EOM / 2 * W[jd, jy] * (g[jh, jd, jy] - g_bar[jh, jd, jy])^2 for jh in JH, jd in JD, jy in JY)
        )

        # Updating CVAR constraint

        if γ < 1
            for jy in JY
                delete(mod, mod.ext[:constraints][:VAR_threshold][jy])
            end
            
            # CVAR constraint
            mod.ext[:constraints][:VAR_threshold] = @constraint(mod, [jy = JY],
            α - profit[jy] <= u[jy] )
        end

    elseif m == "H2turbine"

        VC = mod.ext[:parameters][:VC]
        IC = mod.ext[:parameters][:IC] # overnight investment costs
        η_H2_E = mod.ext[:parameters][:η_H2_E]
        WTP_HCM = 10000 # random, must be chosen

        g = mod.ext[:variables][:g]
        cap = mod.ext[:variables][:cap]
        cap_cm = mod.ext[:variables][:cap_cm]
        capH_cm = mod.ext[:variables][:capH_cm]
        α = mod.ext[:variables][:α]
        u = mod.ext[:variables][:u]

        # Expressions

        #inv_cost = mod.ext[:expressions][:inv_cost] = @expression(mod, IC * cap)

        gH = mod.ext[:expressions][:gH] = @expression(mod, - g / η_H2_E )
        
        cost = mod.ext[:expressions][:cost] = @expression(mod, [jy = JY],
        IC * cap
        + sum(W[jd, jy] * VC * g[jh, jd, jy] for jh in JH, jd in JD)
        - sum(W[jd, jy] * λ_H2[jh, jd, jy] * gH[jh, jd, jy] for jh in JH, jd in JD) )
        
        revenue = mod.ext[:expressions][:revenue] = @expression(mod, [jy = JY],
        σ * λ_CM * cap_cm
        - σH * (WTP_HCM - λ_HCM) * capH_cm   # capH_cm is negative (demand) so there is a - in front of it
        + sum(W[jd, jy] * λ_EOM[jh, jd, jy] * g[jh, jd, jy] for jh in JH, jd in JD) )

        profit = mod.ext[:expressions][:profit] = @expression(mod, [jy = JY],
        revenue[jy] - cost[jy] )

        CVAR = mod.ext[:expressions][:CVAR] = @expression(mod,
        α - ((1/β) * sum(P[jy] * u[jy] for jy in JY)) )

        #tot_revenue = mod.ext[:expressions][:tot_revenue] = @expression(mod,
        #σ * λ_CM * cap_cm + sum(P[jy] * sum(W[jd, jy] * λ_EOM[jh, jd, jy] * g[jh, jd, jy] for jh in JH, jd in JD) for jy in JY))

        # Update objective function 
        mod.ext[:objective] = @objective(mod, Min,
        - γ * sum(P[jy] * (revenue[jy] - cost[jy]) for jy in JY) - (1 - γ) * CVAR
        + sum(ρ_EOM / 2 * W[jd, jy] * (g[jh, jd, jy] - g_bar[jh, jd, jy])^2 for jh in JH, jd in JD, jy in JY)
        + sum(ρ_H2 / 2 * W[jd, jy] * (gH[jh, jd, jy] - gH_bar[jh, jd, jy])^2 for jh in JH, jd in JD, jy in JY)
        + σ * ρ_CM / 2 * (cap_cm - cap_bar)^2
        + σH * ρ_HCM / 2 * (capH_cm - capH_bar)^2)

        # Updating CVAR constraint

        if γ < 1
            for jy in JY
                delete(mod, mod.ext[:constraints][:VAR_threshold][jy])
            end

            # CVAR constraint
            mod.ext[:constraints][:VAR_threshold] = @constraint(mod, [jy = JY],
            α - (revenue[jy] - cost[jy]) <= u[jy] )
        end

    elseif m == "Biomass"

        VC = mod.ext[:parameters][:VC] # variable costs
        IC = mod.ext[:parameters][:IC] # annuity investment costs

        g = mod.ext[:variables][:g]
        cap = mod.ext[:variables][:cap]
        cap_cm = mod.ext[:variables][:cap_cm]
        α = mod.ext[:variables][:α]
        u = mod.ext[:variables][:u]
        
        #inv_cost = mod.ext[:expressions][:inv_cost]
        
        cost = mod.ext[:expressions][:cost] = @expression(mod, [jy = JY],
        IC * cap
        + sum(W[jd, jy] * VC * g[jh, jd, jy] for jh in JH, jd in JD) )

        revenue = mod.ext[:expressions][:revenue] = @expression(mod, [jy = JY],
        σ * λ_CM * cap_cm 
        + sum(W[jd, jy] * λ_EOM[jh, jd, jy] * g[jh, jd, jy] for jh in JH, jd in JD) )

        profit = mod.ext[:expressions][:profit] = @expression(mod, [jy = JY],
        revenue[jy] - cost[jy] )

        CVAR = mod.ext[:expressions][:CVAR] = @expression(mod,
        α - ((1/β) * sum(P[jy] * u[jy] for jy in JY)) )

        #tot_revenue = mod.ext[:expressions][:tot_revenue] = @expression(mod,
        #σ * λ_CM * cap_cm + sum(P[jy] * sum(W[jd, jy] * λ_EOM[jh, jd, jy] * g[jh, jd, jy] for jh in JH, jd in JD) for jy in JY))

        # Update objective function 
        mod.ext[:objective] = @objective(mod, Min,
        - γ * sum(P[jy] * (revenue[jy] - cost[jy]) for jy in JY) - (1 - γ) * CVAR
        + sum(ρ_EOM / 2 * W[jd, jy] * (g[jh, jd, jy] - g_bar[jh, jd, jy])^2 for jh in JH, jd in JD, jy in JY)
        + σ * ρ_CM / 2 * (cap_cm - cap_bar)^2)

        # Updating CVAR constraint

        if γ < 1
            for jy in JY
                delete(mod, mod.ext[:constraints][:VAR_threshold][jy])
            end

            # CVAR constraint
            mod.ext[:constraints][:VAR_threshold] = @constraint(mod, [jy = JY],
            α - (revenue[jy] - cost[jy]) <= u[jy] )
        end


    else
        
        VC = mod.ext[:parameters][:VC]
        IC = mod.ext[:parameters][:IC] # overnight investment costs

        g = mod.ext[:variables][:g]
        cap = mod.ext[:variables][:cap]
        α = mod.ext[:variables][:α]
        u = mod.ext[:variables][:u]
        
        #inv_cost = mod.ext[:expressions][:inv_cost]
        
        cost = mod.ext[:expressions][:cost] = @expression(mod, [jy = JY],
        IC * cap
        + sum(W[jd, jy] * VC * g[jh, jd, jy] for jh in JH, jd in JD) )
        
        revenue = mod.ext[:expressions][:revenue] = @expression(mod, [jy = JY],
        sum(W[jd, jy] * λ_EOM[jh, jd, jy] * g[jh, jd, jy] for jh in JH, jd in JD) )

        profit = mod.ext[:expressions][:profit] = @expression(mod, [jy = JY],
        revenue[jy] - cost[jy] )

        #tot_revenue = mod.ext[:expressions][:tot_revenue] = @expression(mod,
        #sum(P[jy] * revenue[jy] for jy in JY))

        CVAR = mod.ext[:expressions][:CVAR] = @expression(mod,
        α - ((1/β) * sum(P[jy] * u[jy] for jy in JY)) )


        # Update objective function 
        mod.ext[:objective] = @objective(mod, Min,
        - γ * sum(P[jy] * (revenue[jy] - cost[jy]) for jy in JY) - (1 - γ) * CVAR
        + sum(ρ_EOM / 2 * W[jd, jy] * (g[jh, jd, jy] - g_bar[jh, jd, jy])^2 for jh in JH, jd in JD, jy in JY))

        # Updating CVAR constraint

        if γ < 1
            for jy in JY
                delete(mod, mod.ext[:constraints][:VAR_threshold][jy])
            end

            # CVAR constraint
            mod.ext[:constraints][:VAR_threshold] = @constraint(mod, [jy = JY],
            α - (revenue[jy] - cost[jy]) <= u[jy] )
        end

    end

    # solve problem
    optimize!(mod)

    return mod
end