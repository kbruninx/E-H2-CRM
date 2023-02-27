function build_ps_agent!(m::String,mod::Model,EOM::Dict)
    # Extract sets
    JH = mod.ext[:sets][:JH]
    JD = mod.ext[:sets][:JD]

    # Extract parameters
    W = mod.ext[:parameters][:W] # weight of the representative days

    # ADMM algorithm parameters
    λ_H2 = mod.ext[:parameters][:λ_H2] # H2 prices
    gH_bar = mod.ext[:parameters][:gH_bar] # element in ADMM penalty term related to hydrogen market
    ρ_H2 = mod.ext[:parameters][:ρ_H2] # rho-value in ADMM related to H2 market
    λ_EOM = mod.ext[:parameters][:λ_EOM] # EOM prices
    g_bar = mod.ext[:parameters][:g_bar] # element in ADMM penalty term related to EOM
    ρ_EOM = mod.ext[:parameters][:ρ_EOM] # rho-value in ADMM related to EOM auctions

    if m == "H2turbine"
        # Extract time series data
        AF = mod.ext[:timeseries][:AF] # avaiability factors

        # Define parameters and expressions
        
        VC = mod.ext[:parameters][:VC] # variable costs
        IC = mod.ext[:parameters][:IC] # annuity investment costs
        η_H2_E = mod.ext[:parameters][:η_H2_E] # efficiency of hydrogen turbines

        # Create variables
        cap = mod.ext[:variables][:cap] = @variable(mod, lower_bound = 0, base_name = "capacity")
        g = mod.ext[:variables][:g] = @variable(mod, [jh = JH, jd = JD], lower_bound = 0, base_name = "generation")

        # Create expressions
        gH = mod.ext[:expressions][:gH] = @expression(mod,
        - g / η_H2_E )
        tot_cost = mod.ext[:expressions][:tot_cost] = @expression(mod,
        IC * cap + sum(W[jd] * VC * g[jh,jd] for jh in JH, jd in JD) )

        # Objective 
        mod.ext[:objective] = @objective(mod, Min,
        + tot_cost
        - sum(W[jd] * λ_H2[jh, jd] * gH[jh, jd] for jh in JH, jd in JD)
        - sum(W[jd] * λ_EOM[jh, jd] * g[jh, jd] for jh in JH, jd in JD)
        + sum(ρ_EOM / 2 * W[jd] * (g[jh, jd] - g_bar[jh, jd])^2 for jh in JH, jd in JD)
        + sum(ρ_H2 / 2 * W[jd] * (gH[jh,jd] - gH_bar[jh,jd])^2 for jh in JH, jd in JD))

        # Capacity constraint
        mod.ext[:constraints][:cap_limit] = @constraint(mod, [jh = JH, jd = JD],
        g[jh, jd] <= AF[jh, jd] * cap / 1000 )  # scaling factor needed to go from GW -> TWh

    elseif m == "Edemand"
        # Define parameters and expressions
        WTP = mod.ext[:parameters][:WTP]  # willingness to pay ("price cap") of consumers
        ela = EOM["elasticity"] # section of price-elasticity [MWh]

        # Create variables
        g_VOLL = mod.ext[:variables][:g_VOLL] = @variable(mod, [jh = JH, jd = JD], upper_bound = 0, base_name = "generation")
        g_ela = mod.ext[:variables][:g_ela] = @variable(mod, [jh = JH, jd = JD], upper_bound = 0, base_name = "generation")
        # g for demand agent is defined as negative

        # Create expressions
        g =  mod.ext[:expressions][:g] = @expression(mod, g_VOLL + g_ela)

        # Objective 
        mod.ext[:objective] = @objective(mod, Min,   # minimize a negative quantity (maximize its absolute value)
        sum(W[jd] * ((WTP * g[jh,jd] + (g_ela[jh,jd])^2 * WTP / (2*ela[jh,jd])) - λ_EOM[jh, jd] * g[jh,jd]) for jh in JH, jd in JD)
        + sum(ρ_EOM / 2 * W[jd] * (g[jh, jd] - g_bar[jh, jd])^2 for jh in JH, jd in JD)
        ) 
        # g_ela is defined negative, therefore the function here must be adapted since there is g_ela^2 : put + instead of - in front of g_ela^2
        # this definition allows to firstly maximize g_voll and than g_ela, because g_voll gives more 

        # Constraints
        mod.ext[:constraints][:elastic_demand_limit] = @constraint(mod, [jh = JH, jd = JD],
        0.2 * EOM["D"][jh, jd] + g_ela[jh,jd] >= 0)   # again, g is negative
        
        mod.ext[:constraints][:demand_limit] = @constraint(mod, [jh = JH, jd = JD], 
        0.8 * EOM["D"][jh, jd] + g_VOLL[jh, jd] >= 0)

    else
        # Extract parameters
        VC = mod.ext[:parameters][:VC] # variable costs
        IC = mod.ext[:parameters][:IC] # annuity investment costs

        # Create variables
        cap = mod.ext[:variables][:cap] = @variable(mod, lower_bound = 0, base_name = "capacity")
        g = mod.ext[:variables][:g] = @variable(mod, [jh = JH, jd = JD], lower_bound = 0, base_name = "generation")
        tot_cost = mod.ext[:expressions][:tot_cost] = @expression(mod,
        IC * cap + sum(W[jd] * VC * g[jh, jd] for jh in JH, jd in JD))

        # Extract time series data
        AF = mod.ext[:timeseries][:AF] # avaiability factors

        # Objective 
        mod.ext[:objective] = @objective(mod, Min,
        + tot_cost
        - sum(W[jd] * λ_EOM[jh, jd] * g[jh, jd] for jh in JH, jd in JD)
        + sum(ρ_EOM / 2 * W[jd] * (g[jh, jd] - g_bar[jh, jd])^2 for jh in JH, jd in JD))

        # Capacity constraint
        mod.ext[:constraints][:cap_limit] = @constraint(mod, [jh = JH, jd = JD],
        g[jh, jd] <= AF[jh, jd] * cap / 1000 )  # scaling factor needed to go from GW -> TWh
    end   

    return mod
end