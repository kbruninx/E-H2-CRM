function build_h2s_agent!(m::String,mod::Model)
    # Extract sets
    JH = mod.ext[:sets][:JH]
    JD = mod.ext[:sets][:JD]
       
    # Extract common parameters
    W = mod.ext[:parameters][:W] # weight of the representative days
    IC = mod.ext[:parameters][:IC] # annuity investment costs

    # ADMM algorithm parameters
    λ_EOM = mod.ext[:parameters][:λ_EOM] # EOM prices
    g_bar = mod.ext[:parameters][:g_bar] # element in ADMM penalty term related to EOM
    ρ_EOM = mod.ext[:parameters][:ρ_EOM] # rho-value in ADMM related to EUA auctions
    λ_H2 = mod.ext[:parameters][:λ_H2] # H2 prices
    gH_bar = mod.ext[:parameters][:gH_bar] # element in ADMM penalty term related to hydrogen market
    ρ_H2 = mod.ext[:parameters][:ρ_H2] # rho-value in ADMM related to H2 market
    

    if m == "electrolysis" 

        # Extract parameters
        η_E_H2 = mod.ext[:parameters][:η_E_H2] # efficiency electrolysis

        # Decision variables
        capH = mod.ext[:variables][:capH] = @variable(mod, lower_bound=0, base_name="capacity")
        g = mod.ext[:variables][:g] = @variable(mod, [jh=JH,jd=JD], upper_bound=0, base_name="demand_electricity_hydrogen") # note this is defined as a negative number, consumption

        # Create affine expressions  
        inv_cost = mod.ext[:expressions][:inv_cost] = @expression(mod, IC*capH)
        gH = mod.ext[:expressions][:gH] = @expression(mod, -η_E_H2 * g) # [TWh]

        # Definition of the objective function
        mod.ext[:objective] = @objective(mod, Min,
        + inv_cost # [MEUR]
        - sum(W[jd]*(λ_EOM[jh,jd])*g[jh,jd] for jh in JH, jd in JD) # [MEUR] (λ in [MEUR/TWh]=[EUR/MWh])
        - sum(W[jd]*λ_H2[jh,jd]*gH[jh,jd] for jh in JH, jd in JD)
        + sum(ρ_EOM/2*W[jd]*(g[jh,jd] - g_bar[jh,jd])^2 for jh in JH, jd in JD)  
        + sum(ρ_H2/2*W[jd]*(gH[jh,jd] - gH_bar[jh,jd])^2 for jh in JH, jd in JD) 
        )

        # Constraints

        # Electricity consumption
        mod.ext[:constraints][:elec_consumption] = @constraint(mod, [jh=JH,jd=JD], 
        -η_E_H2*g[jh,jd] <= capH/1000  # [TWh]
        )

    elseif m == "H2storage"
        
        # Extract parameters
        VC = mod.ext[:parameters][:VC] # variable cost storage
        η_ch = mod.ext[:parameters][:η_ch] # charging efficiency 
        η_dh = mod.ext[:parameters][:η_dh] # discharging efficiency

        # Decision variables
        capH = mod.ext[:variables][:capH] = @variable(mod, lower_bound=0, base_name="capacity")
        volH = mod.ext[:variables][:volH] = @variable(mod, lower_bound=0, base_name="volume")
        chH = mod.ext[:variables][:chH] = @variable(mod, [jh=JH,jd=JD], lower_bound=0, base_name="hydrogen_charging")
        dhH = mod.ext[:variables][:dhH] = @variable(mod, [jh=JH,jd=JD], lower_bound=0, base_name="hydrogen_discharging")
        SOC = mod.ext[:variables][:SOC] = @variable(mod, [jh=JH,jd=JD], lower_bound=0, base_name="state_of_charge")


        # Create affine expressions  
        inv_cost = mod.ext[:expressions][:inv_cost] = @expression(mod, IC*capH) # add a term for volume in overview_data
        # SOC = mod.ext[:expressions][:SOC] = @expression(mod, SOC + η_ch*chH - dhH/η_dh) # state of charge update
        gH = mod.ext[:expressions][:gH]  = @expression(mod, dhH-chH);

        # Definition of the objective function
        mod.ext[:objective] = @objective(mod, Min,
        + inv_cost # [MEUR]
        - sum(W[jd]*(λ_H2[jh,jd])*dhH[jh,jd] for jh in JH, jd in JD) # [MEUR]
        + sum(W[jd]*λ_H2[jh,jd]*chH[jh,jd] for jh in JH, jd in JD)
        + sum(ρ_H2/2*W[jd]*(gH[jh,jd] - gH_bar[jh,jd])^2 for jh in JH, jd in JD) 
        )

        # Constraints
        mod.ext[:constraints][:dischargin_lim] = @constraint(mod, [jh=JH,jd=JD], 
        dhH[jh,jd] <= capH/1000 )  # [TWh]
        
        mod.ext[:constraints][:charging_lim] = @constraint(mod, [jh=JH,jd=JD], 
        chH[jh,jd] <= capH/1000 ) # [TWh]

        mod.ext[:constraints][:storage_balance] = @constraint(mod,
        sum(W[jd]*(η_ch*chH[jh,jd] - dhH[jh,jd]/η_dh) for jh in JH, jd in JD) == 0)

        mod.ext[:constraints][:SOC_update] = @constraint(mod, [jh=first(JH),jd=JD[2:1:end]],    # divided into two constraints to consider the first hour of the day that has to refer to the last one of the previous day
        SOC[jh,jd] == SOC[last(JH),jd-1] + η_ch*chH[jh,jd] - dhH[jh,jd]/η_dh )

        mod.ext[:constraints][:SOC_update] = @constraint(mod, [jh=JH[2:1:end],jd=JD], 
        SOC[jh, jd] == SOC[jh-1, jd] + η_ch * chH[jh, jd] - dhH[jh, jd] / η_dh )    
    
        mod.ext[:constraints][:SOC_limit] = @constraint(mod, [jh = JH, jd = JD],
        SOC[jh, jd] <= volH)
       
        mod.ext[:constraints][:SOC_0] = @constraint(mod,
        SOC[1,1] == 0)

    end
    
    return mod

end