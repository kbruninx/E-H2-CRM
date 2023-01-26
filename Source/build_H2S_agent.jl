function build_h2s_agent!(mod::Model)
    # Extract sets
    JH = mod.ext[:sets][:JH]
    JD = mod.ext[:sets][:JD]
       
    # Extract parameters
    W = mod.ext[:parameters][:W] # weight of the representative days
    IC = mod.ext[:parameters][:IC] # annuity investment costs
    η_E_H2 = mod.ext[:parameters][:η_E_H2] # efficiency E->H2

    # ADMM algorithm parameters
    λ_EOM = mod.ext[:parameters][:λ_EOM] # EOM prices
    g_bar = mod.ext[:parameters][:g_bar] # element in ADMM penalty term related to EOM
    ρ_EOM = mod.ext[:parameters][:ρ_EOM] # rho-value in ADMM related to EUA auctions
    λ_H2 = mod.ext[:parameters][:λ_H2] # H2 prices
    gH_bar = mod.ext[:parameters][:gH_bar] # element in ADMM penalty term related to hydrogen market
    ρ_H2 = mod.ext[:parameters][:ρ_H2] # rho-value in ADMM related to H2 market
    
    # Decision variables
    capH = mod.ext[:variables][:capH] = @variable(mod, lower_bound=0, base_name="capacity")
    gH = mod.ext[:variables][:gH] = @variable(mod, [jh=JH,jd=JD], lower_bound=0, base_name="generation_hydrogen")
    g = mod.ext[:variables][:g] = @variable(mod, [jh=JH,jd=JD], upper_bound=0, base_name="demand_electricity_hydrogen") # note this is defined as a negative number, consumption

    # Create affine expressions  
    mod.ext[:expressions][:tot_cost] = @expression(mod, 
        IC*capH 
    )

    # Definition of the objective function
    mod.ext[:objective] = @objective(mod, Min,
    + IC*capH # [MEUR]
    - sum(W[jd]*(λ_EOM[jh,jd])*g[jh,jd] for jh in JH, jd in JD) # [MEUR]
    - sum(W[jd]*λ_H2[jh,jd]*gH[jh,jd] for jh in JH, jd in JD)
    + sum(ρ_EOM/2*W[jd]*(g[jh,jd] - g_bar[jh,jd])^2 for jh in JH, jd in JD)  
    + sum(ρ_H2/2*W[jd]*(gH[jh,jd] - gH_bar[jh,jd])^2 for jh in JH, jd in JD) 
    )
    
    # Constraints
    mod.ext[:constraints][:gen_limit_capacity] = @constraint(mod, [jh=JH,jd=JD],
        gH[jh,jd] <=  capH/1000 #  [TWh]        
    )

    # Electricity consumption
    mod.ext[:constraints][:elec_consumption] = @constraint(mod, [jh=JH,jd=JD],
        -η_E_H2*g[jh,jd] <= capH/1000  # [TWh]
    )    

    # H2 production
    mod.ext[:constraints][:gen_limit_energy_sources] = @constraint(mod, [jh=JH,jd=JD],
        gH[jh,jd] <= -η_E_H2*g[jh,jd]  # [TWh]
    )
 
    return mod

end