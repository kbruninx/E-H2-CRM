function build_ps_agent!(m::String,mod::Model,EOM::Dict)
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
    λ_H2 = mod.ext[:parameters][:λ_H2] # H2 prices
    gH_bar = mod.ext[:parameters][:gH_bar] # element in ADMM penalty term related to hydrogen market
    ρ_H2 = mod.ext[:parameters][:ρ_H2] # rho-value in ADMM related to H2 market
    λ_EOM = mod.ext[:parameters][:λ_EOM] # EOM prices
    g_bar = mod.ext[:parameters][:g_bar] # element in ADMM penalty term related to EOM
    ρ_EOM = mod.ext[:parameters][:ρ_EOM] # rho-value in ADMM related to EOM auctions
    λ_CM = mod.ext[:parameters][:λ_CM] # capacity market prices
    cap_bar = mod.ext[:parameters][:cap_bar]  # element in ADMM penalty term related to CM
    ρ_CM = mod.ext[:parameters][:ρ_CM]  # rho-value in ADMM related to CM auctions
    λ_HCM = mod.ext[:parameters][:λ_HCM] # hydrogen capacity market prices
    capH_bar = mod.ext[:parameters][:capH_bar]  # element in ADMM penalty term related to hydrogen CM
    ρ_HCM = mod.ext[:parameters][:ρ_HCM]  # rho-value in ADMM related to hydrogen CM auctions


    if m == "Edemand"
        # Define parameters and expressions
        WTP = mod.ext[:parameters][:WTP]  # willingness to pay ("price cap") of consumers
        ela = EOM["elasticity"] # section of price-elasticity [10^2 GWh]

        # Create variables
        g_VOLL = mod.ext[:variables][:g_VOLL] = @variable(mod, [jh = JH, jd = JD, jy = JY], lower_bound = 0, base_name = "demand")
        g_ela = mod.ext[:variables][:g_ela] = @variable(mod, [jh = JH, jd = JD, jy = JY], lower_bound = 0, base_name = "elastic demand")
        α = mod.ext[:variables][:α] = @variable(mod, base_name = "VAR")
        u = mod.ext[:variables][:u] = @variable(mod, [jy = JY], lower_bound = 0, base_name = "tail profit difference")  # profit difference of worst-case tail scenarios with respect to the VAR

        # Warm start
        #=
        for jy in JY
            start_value = Matrix(CSV.read(joinpath(home_dir, "Input", "RA10_elCM", string("Elastic_demand_el_", jy, ".csv")), delim=";", DataFrame))
            for jh in JH, jd in JD, jy in JY
                set_start_value(g_ela[jh, jd, jy], start_value[jh, jd])
            end
        end
        =#

        # Create expressions
        g_positive = mod.ext[:expressions][:g_positive] = @expression(mod, g_VOLL + g_ela)
        g =  mod.ext[:expressions][:g] = @expression(mod, - g_positive)

        profit = mod.ext[:expressions][:profit] = @expression(mod, [jy = JY],      # the profit is formulated negative to make the problem convex
        sum(W[jd, jy] * (WTP * g_positive[jh, jd, jy] - (g_ela[jh, jd, jy])^2 * WTP / (2*ela[jh, jd, jy]) 
        - λ_EOM[jh, jd, jy] * g_positive[jh, jd, jy]) for jh in JH, jd in JD)
        - σ * λ_CM * CM["D"] )

        #revenue = mod.ext[:expressions][:revenue] = @expression(mod, [jy = JY],
        #sum(W[jd, jy] * (WTP * g_positive[jh, jd, jy] - (g_ela[jh, jd, jy])^2 * WTP / (2*ela[jh, jd, jy])) for jh in JH, jd in JD) )

        cost = mod.ext[:expressions][:cost] = @expression(mod, [jy = JY],
        sum(W[jd, jy] * λ_EOM[jh, jd, jy] * g_positive[jh, jd, jy] for jh in JH, jd in JD)
        + σ * λ_CM * CM["D"] )

        CVAR = mod.ext[:expressions][:CVAR] = @expression(mod,
        α - ((1/β) * sum(P[jy] * u[jy] for jy in JY)) )

        # Objective 
        mod.ext[:objective] = @objective(mod, Min,   # minimize a negative quantity (maximize its absolute value)
        - γ * sum(P[jy] * profit[jy] for jy in JY) - (1 - γ) * CVAR                   # profit is already intended as negative (as cost-revenue), so no need for - in front of it
        + sum(ρ_EOM / 2 * W[jd, jy] * (g[jh, jd, jy] - g_bar[jh, jd, jy])^2 for jh in JH, jd in JD, jy in JY)
        )
        
        # Constraints
        mod.ext[:constraints][:elastic_demand_limit] = @constraint(mod, [jh = JH, jd = JD, jy = JY],
        g_ela[jh, jd, jy] <= ela[jh, jd, jy]) 
        
        mod.ext[:constraints][:demand_limit] = @constraint(mod, [jh = JH, jd = JD, jy = JY], 
        g_VOLL[jh, jd, jy] <= 0.8 * EOM["D"][jh, jd, jy] )

        if γ < 1
            # CVAR constraint
            mod.ext[:constraints][:VAR_threshold] = @constraint(mod, [jy = JY],
            α - profit[jy] <= u[jy]) 
        end        


    elseif m == "H2turbine"
        # Extract time series data
        AF = mod.ext[:timeseries][:AF] # avaiability factors

        # Define parameters and expressions
        
        VC = mod.ext[:parameters][:VC] # variable costs
        IC = mod.ext[:parameters][:IC] # annuity investment costs
        η_H2_E = mod.ext[:parameters][:η_H2_E] # efficiency of hydrogen turbines

        # Create variables
        cap = mod.ext[:variables][:cap] = @variable(mod, lower_bound = 0, base_name = "capacity")
        capH_cm = mod.ext[:variables][:capH_cm] = @variable(mod, upper_bound = 0, base_name = "capacity demand HCM")
        g = mod.ext[:variables][:g] = @variable(mod, [jh = JH, jd = JD, jy = JY], lower_bound = 0, base_name = "generation")
        cap_cm = mod.ext[:variables][:cap_cm] = @variable(mod, lower_bound = 0, base_name = "capacity offered CM")
        α = mod.ext[:variables][:α] = @variable(mod, base_name = "VAR")
        u = mod.ext[:variables][:u] = @variable(mod, [jy = JY], lower_bound = 0, base_name = "tail profit difference")  # profit difference of worst-case tail scenarios with respect to the VAR


        # Warm start
        #=
        for jy in JY
            start_value = Matrix(CSV.read(joinpath(home_dir, "Input", "RA10_elCM", string("GenerationH2turbine_", jy, ".csv")), delim=";", DataFrame))
            for jh in JH, jd in JD, jy in JY
                set_start_value(g[jh, jd, jy], start_value[jh, jd])
            end
        end
        set_start_value(cap, (CSV.read(joinpath(home_dir, "Input", "RA10_elCM", "capacity_ps.csv"), delim=";", DataFrame))[!,:CAP_H2turbine][1])
        =#

        # Create expressions
        gH = mod.ext[:expressions][:gH] = @expression(mod, - g / η_H2_E )
        
        # inv_cost = mod.ext[:expressions][:inv_cost] = @expression(mod, IC * cap)
        
        cost = mod.ext[:expressions][:cost] = @expression(mod, [jy = JY],
        IC * cap
        + sum(W[jd, jy] * VC * g[jh, jd, jy] for jh in JH, jd in JD)
        - sum(W[jd, jy] * λ_H2[jh, jd, jy] * gH[jh, jd, jy] for jh in JH, jd in JD)
        - σH * λ_HCM * capH_cm)   # capH_cm is negative (demand) so there is a - in front of it 
        
        revenue = mod.ext[:expressions][:revenue] = @expression(mod, [jy = JY],
        σ * λ_CM * cap_cm
        + sum(W[jd, jy] * λ_EOM[jh, jd, jy] * g[jh, jd, jy] for jh in JH, jd in JD) )

        profit = mod.ext[:expressions][:profit] = @expression(mod, [jy = JY],
        revenue[jy] - cost[jy]  )

        CVAR = mod.ext[:expressions][:CVAR] = @expression(mod,
        α - ((1/β) * sum(P[jy] * u[jy] for jy in JY)) )

        # Objective 
        mod.ext[:objective] = @objective(mod, Min,
        - γ * sum(P[jy] * (revenue[jy] - cost[jy]) for jy in JY) - (1 - γ) * CVAR
        + sum(ρ_EOM / 2 * W[jd, jy] * (g[jh, jd, jy] - g_bar[jh, jd, jy])^2 for jh in JH, jd in JD, jy in JY)
        + sum(ρ_H2 / 2 * W[jd, jy] * (gH[jh, jd, jy] - gH_bar[jh, jd, jy])^2 for jh in JH, jd in JD, jy in JY)
        + σ * ρ_CM / 2 * (cap_cm - cap_bar)^2
        + σH * ρ_HCM / 2 * (capH_cm - capH_bar)^2)

        # Capacity constraint
        mod.ext[:constraints][:cap_limit] = @constraint(mod, [jh = JH, jd = JD, jy = JY],
        g[jh, jd, jy] <= AF[jh, jd, jy] * cap / 100)

        if σ == 1
            mod.ext[:constraints][:CM] = @constraint(mod,
            cap_cm <= cap )
        end

        if σH == 1  
            mod.ext[:constraints][:HCM] = @constraint(mod,
            - capH_cm == cap_cm / η_H2_E)    #### H2 turbine must offer firm capacity
        end

        if γ < 1
            # CVAR constraint
            mod.ext[:constraints][:VAR_threshold] = @constraint(mod, [jy = JY],
            α - (revenue[jy] - cost[jy]) <= u[jy] )
        end
        

    elseif m == "Biomass"

        # Extract parameters
        VC = mod.ext[:parameters][:VC] # variable costs
        IC = mod.ext[:parameters][:IC] # annuity investment costs

        # Create variables
        cap = mod.ext[:variables][:cap] = @variable(mod, lower_bound = 0, base_name = "capacity")
        g = mod.ext[:variables][:g] = @variable(mod, [jh = JH, jd = JD, jy = JY], lower_bound = 0, base_name = "generation")
        cap_cm = mod.ext[:variables][:cap_cm] = @variable(mod, lower_bound = 0, base_name = "capacity offered")
        α = mod.ext[:variables][:α] = @variable(mod, base_name = "VAR")
        u = mod.ext[:variables][:u] = @variable(mod, [jy = JY], lower_bound = 0, base_name = "tail profit difference")  # profit difference of worst-case tail scenarios with respect to the VAR
        
        # Warm start 
        #=
        for jy in JY
            start_value = Matrix(CSV.read(joinpath(home_dir, "Input", "RA10_elCM", string("GenerationBiomass_", jy, ".csv")), delim=";", DataFrame))
            for jh in JH, jd in JD, jy in JY
                set_start_value(g[jh, jd, jy], start_value[jh, jd])
            end
        end
        set_start_value(cap, (CSV.read(joinpath(home_dir, "Input", "RA10_elCM", "capacity_ps.csv"), delim=";", DataFrame))[!, :CAP_Biomass][1])
        =#

        # Create expressions
        #inv_cost = mod.ext[:expressions][:inv_cost] = @expression(mod, IC * cap )

        revenue = mod.ext[:expressions][:revenue] = @expression(mod, [jy = JY],
        σ * λ_CM * cap_cm 
        + sum(W[jd, jy] * λ_EOM[jh, jd, jy] * g[jh, jd, jy] for jh in JH, jd in JD) )

        cost = mod.ext[:expressions][:cost] = @expression(mod, [jy = JY],
        IC * cap
        + sum(W[jd, jy] * VC * g[jh, jd, jy] for jh in JH, jd in JD) )

        CVAR = mod.ext[:expressions][:CVAR] = @expression(mod,
        α - ((1/β) * sum(P[jy] * u[jy] for jy in JY)) )

        profit = mod.ext[:expressions][:profit] = @expression(mod, [jy = JY],
        revenue[jy] - cost[jy] )

        # Extract time series data
        AF = mod.ext[:timeseries][:AF] # avaiability factors

        # Objective 
        mod.ext[:objective] = @objective(mod, Min,
        - γ * sum(P[jy] * (revenue[jy] - cost[jy]) for jy in JY) - (1 - γ) * CVAR
        + sum(ρ_EOM / 2 * W[jd, jy] * (g[jh, jd, jy] - g_bar[jh, jd, jy])^2 for jh in JH, jd in JD, jy in JY)
        + σ * ρ_CM / 2 * (cap_cm - cap_bar)^2 )

        # Capacity constraint
        mod.ext[:constraints][:cap_limit] = @constraint(mod, [jh = JH, jd = JD, jy = JY],
        g[jh, jd, jy] <= AF[jh, jd, jy] * cap / 100)

        if σ == 1
            mod.ext[:constraints][:CM] = @constraint(mod,
            cap_cm <= cap ) 
        end

        if γ < 1
            # CVAR constraint
            mod.ext[:constraints][:VAR_threshold] = @constraint(mod, [jy = JY],
            α - (revenue[jy] - cost[jy]) <= u[jy] )
        end


    else

        # Extract parameters
        VC = mod.ext[:parameters][:VC] # variable costs
        IC = mod.ext[:parameters][:IC] # annuity investment costs

        # Create variables
        cap = mod.ext[:variables][:cap] = @variable(mod, lower_bound = 0, base_name = "capacity")
        g = mod.ext[:variables][:g] = @variable(mod, [jh = JH, jd = JD, jy = JY], lower_bound = 0, base_name = "generation")
        α = mod.ext[:variables][:α] = @variable(mod, base_name = "VAR")
        u = mod.ext[:variables][:u] = @variable(mod, [jy = JY], lower_bound = 0, base_name = "tail profit difference")  # profit difference of worst-case tail scenarios with respect to the VAR
        
        # Warm start
        #for jy in JY
        #    start_value = Matrix(CSV.read(joinpath(home_dir, "Input", "RA10_elCM", string("Generation", m, "_", jy, ".csv")), delim=";", DataFrame))
        #    for jh in JH, jd in JD, jy in JY
        #        set_start_value(g[jh, jd, jy], start_value[jh, jd])
        #    end
        #end

        # Create expressions
        # inv_cost = mod.ext[:expressions][:inv_cost] = @expression(mod, IC * cap )

        revenue = mod.ext[:expressions][:revenue] = @expression(mod, [jy = JY],
        sum(W[jd, jy] * λ_EOM[jh, jd, jy] * g[jh, jd, jy] for jh in JH, jd in JD) )

        cost = mod.ext[:expressions][:cost] = @expression(mod, [jy = JY],
        IC * cap
        + sum(W[jd, jy] * VC * g[jh, jd, jy] for jh in JH, jd in JD) )

        profit = mod.ext[:expressions][:profit] = @expression(mod, [jy = JY],
        revenue[jy] - cost[jy] )

        CVAR = mod.ext[:expressions][:CVAR] = @expression(mod,
        α - ((1/β) * sum(P[jy] * u[jy] for jy in JY)) )

        # Extract time series data
        AF = mod.ext[:timeseries][:AF] # avaiability factors

        # Objective 
        mod.ext[:objective] = @objective(mod, Min,
        - γ * sum(P[jy] * (revenue[jy] - cost[jy]) for jy in JY) - (1 - γ) * CVAR
        + sum(ρ_EOM / 2 * W[jd, jy] * (g[jh, jd, jy] - g_bar[jh, jd, jy])^2 for jh in JH, jd in JD, jy in JY) )

        # Capacity constraint
        mod.ext[:constraints][:cap_limit] = @constraint(mod, [jh = JH, jd = JD, jy = JY],
        g[jh, jd, jy] <= AF[jh, jd, jy] * cap / 100 )

        #if m == "Solar"
        #    mod.ext[:constraints][:cap_solar] = @constraint(mod,
        #        cap == 71.137)
        #elseif m == "WindOnshore"
        #    mod.ext[:constraints][:cap_Wind] = @constraint(mod,
        #        cap == 69.414)
        #end

        if γ < 1
            # CVAR constraint
            mod.ext[:constraints][:VAR_threshold] = @constraint(mod, [jy = JY],
            α - (revenue[jy] - cost[jy]) <= u[jy] )
        end

    end   

    return mod
end