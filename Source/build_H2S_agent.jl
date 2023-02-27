function build_h2s_agent!(m::String,mod::Model,H2::Dict)
    # Extract sets
    JH = mod.ext[:sets][:JH]
    JD = mod.ext[:sets][:JD]
       
    # Extract common parameters
    W = mod.ext[:parameters][:W] # weight of the representative days

    # ADMM algorithm parameters
    λ_EOM = mod.ext[:parameters][:λ_EOM] # EOM prices
    g_bar = mod.ext[:parameters][:g_bar] # element in ADMM penalty term related to EOM
    ρ_EOM = mod.ext[:parameters][:ρ_EOM] # rho-value in ADMM related to EUA auctions
    λ_H2 = mod.ext[:parameters][:λ_H2] # H2 prices
    gH_bar = mod.ext[:parameters][:gH_bar] # element in ADMM penalty term related to hydrogen market
    ρ_H2 = mod.ext[:parameters][:ρ_H2] # rho-value in ADMM related to H2 market
    

    if m == "H2demand"
        # Extract parameters
        WTP = mod.ext[:parameters][:WTP]  # willingness to pay ("price cap") of consumers
        ela_H2 = H2["elasticity"] # section of price-elasticity [MWh]

        # Create variables
        gH_VOLL = mod.ext[:variables][:gH_VOLL] = @variable(mod, [jh = JH, jd = JD], upper_bound = 0, base_name = "H2_inelastic_demand")
        gH_ela = mod.ext[:variables][:gH_ela] = @variable(mod, [jh = JH, jd = JD], upper_bound = 0, base_name = "H2_elastic_demand")
        # gH for demand agent is defined as negative

        # Create expressions
        gH =  mod.ext[:expressions][:gH] = @expression(mod, gH_VOLL + gH_ela)

        # Objective 
        mod.ext[:objective] = @objective(mod, Min,   # minimize a negative quantity (maximize its absolute value)
        sum(W[jd] * ((WTP * gH[jh,jd] + (gH_ela[jh,jd])^2 * WTP / (2*ela_H2)) - λ_H2[jh, jd] * gH[jh,jd]) for jh in JH, jd in JD)
        + sum(ρ_H2 / 2 * W[jd] * (gH[jh, jd] - gH_bar[jh, jd])^2 for jh in JH, jd in JD)
        )  # N.B. gH is negative here
        # gH_ela is defined negative, therefore the function here must be adapted since there is gH_ela^2 : put + instead of - in front of gH_ela^2

        # Constraints
        mod.ext[:constraints][:elastic_demand_limit] = @constraint(mod, [jh = JH, jd = JD],
        0.2 * H2["D"][jh, jd] + gH_ela[jh,jd] >= 0)   # again, g is negative
        
        mod.ext[:constraints][:demand_limit] = @constraint(mod, [jh = JH, jd = JD], 
        0.8 * H2["D"][jh, jd] + gH_VOLL[jh, jd] >= 0)

    elseif m == "electrolysis" 

        # Extract parameters
        η_E_H2 = mod.ext[:parameters][:η_E_H2] # efficiency electrolysis
        IC = mod.ext[:parameters][:IC] # annuity investment costs

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
        JD_A = mod.ext[:sets][:JD_A] # sets of all days used for storage SOC
        index_repr =  mod.ext[:sets][:index_repr] # index of representative days
        V = mod.ext[:sets][:order_matrix]
        VC = mod.ext[:parameters][:VC] # variable cost storage
        η_ch = mod.ext[:parameters][:η_ch] # charging efficiency 
        η_dh = mod.ext[:parameters][:η_dh] # discharging efficiency
        IC_cap = mod.ext[:parameters][:IC_cap] # annuity investment cost for capacity
        IC_vol = mod.ext[:parameters][:IC_vol] # annuity investment cost for volume

        # Decision variables
        capH = mod.ext[:variables][:capH] = @variable(mod, lower_bound=0, base_name="capacity") # GW
        volH = mod.ext[:variables][:volH] = @variable(mod, lower_bound=0, base_name="volume") # TWh
        chH = mod.ext[:variables][:chH] = @variable(mod, [jh=JH,jd=JD], lower_bound=0, base_name="hydrogen_charging")
        dhH = mod.ext[:variables][:dhH] = @variable(mod, [jh=JH,jd=JD], lower_bound=0, base_name="hydrogen_discharging")
        delta_e = mod.ext[:expressions][:delta_e] = @variable(mod, [jd=JD])
        SOC = mod.ext[:variables][:SOC] = @variable(mod, [jh=JH,jd=JD], lower_bound=0, base_name="state_of_charge")
        SOC_0 = mod.ext[:variables][:SOC_0] = @variable(mod, [jd=JD], lower_bound=0, base_name="initial_state_of_charge_repr_day")    # state of charge at the start of a period
        SOC_AD_0 = mod.ext[:variables][:SOC_AD_0] = @variable(mod, [jd_a=JD_A], lower_bound=0, base_name="initial_state_of_charge_all_days")
        max_deviation_SOC = mod.ext[:variables][:max_deviation_SOC] = @variable(mod, [jd=JD], lower_bound=0, base_name="max_deviaion_SOC")
        min_deviation_SOC = mod.ext[:variables][:min_deviation_SOC] = @variable(mod, [jd=JD], lower_bound=0, base_name="min_deviaion_SOC")    # refer to a negative deviation
        max_deviation_SOC_AD = mod.ext[:variables][:max_deviation_SOC_AD] = @variable(mod, [jd_a=JD_A], lower_bound=0)
        min_deviation_SOC_AD = mod.ext[:variables][:min_deviation_SOC_AD] = @variable(mod, [jd_a=JD_A], lower_bound=0)

        # Create affine expressions  
        inv_cost = mod.ext[:expressions][:inv_cost] = @expression(mod, IC_cap*capH + IC_vol*volH)
        gH = mod.ext[:expressions][:gH]  = @expression(mod, dhH-chH);

        # Definition of the objective function
        mod.ext[:objective] = @objective(mod, Min,
        + inv_cost # [MEUR]
        - sum(W[jd]*(λ_H2[jh,jd])*dhH[jh,jd] for jh in JH, jd in JD) # [MEUR]
        + sum(W[jd]*λ_H2[jh,jd]*chH[jh,jd] for jh in JH, jd in JD)
        + sum(ρ_H2/2*W[jd]*(gH[jh,jd] - gH_bar[jh,jd])^2 for jh in JH, jd in JD) 
        )

        # Constraints

        # Inspired by Gonzato, Bruninx and Delarue - Long term storage in generation expansion planning models with a reduced temporal scope
        #### take in input V (ordering_variable), check dimensions of V in the constraints

        # Linking SOC_0 and SOC_AD_0 
        mod.ext[:constraints][:D2AD] = @constraint(mod, [jd=JD, jd_a=JD_A],
        SOC_0[jd] == SOC_AD_0[Int(index_repr[jd])] )

        # Define energy change for each representative day - (33)
        mod.ext[:constraints][:delta_e_repr] = @constraint(mod, [jh=JH, jd=JD], 
        delta_e[jd] == sum((η_ch*chH[jh,jd] - dhH[jh,jd]/η_dh) for jh in JH) )

        # Base state of charge definition (for all days) - (34) - shouldn't they be of the previous jd_a ?
        mod.ext[:constraints][:SOC_AD_0] = @constraint(mod, [jd=JD, jd_a=JD_A[2:end]],
        SOC_AD_0[jd_a] == SOC_AD_0[jd_a-1] + sum(V[jd_a,jd]*delta_e[jd] for jd in JD) )

        # Max positive and negative deviations from the base state of charge for representiative periods  - (35),(36)
        mod.ext[:constraints][:max_deviation_repr] = @constraint(mod, [jh=JH, jd=JD],
        max_deviation_SOC[jd] >= SOC[jh,jd] - SOC_0[jd] )

        mod.ext[:constraints][:min_deviation_repr] = @constraint(mod, [jh=JH, jd=JD],
        min_deviation_SOC[jd] >=  SOC_0[jd] - SOC[jh,jd] )

        # Max positive and negative deviations from the base state of charge for non-representiative periods - (37),(38)
        mod.ext[:constraints][:max_deviation_all] = @constraint(mod, [jd=JD, jd_a=JD_A],
        max_deviation_SOC_AD[jd_a] == sum(V[jd_a,jd] * max_deviation_SOC[jd] for jd in JD) )

        mod.ext[:constraints][:min_deviation_all] = @constraint(mod, [jd=JD, jd_a=JD_A],
        min_deviation_SOC_AD[jd_a] == sum(V[jd_a,jd] * min_deviation_SOC[jd] for jd in JD) )

        # SOC limits for all periods - (39),(40)
        mod.ext[:constraints][:SOC_upper_limit_all] = @constraint(mod, [jd_a=JD_A],
        volH >= SOC_AD_0[jd_a] + max_deviation_SOC_AD[jd_a] )

        mod.ext[:constraints][:SOC_lower_limit_all] = @constraint(mod, [jd_a=JD_A],
        SOC_AD_0[jd_a] - min_deviation_SOC_AD[jd_a] >= 0 )

        # Cyclic constraint - (41) but adpated, the one in the paper does not make sense
        mod.ext[:constraints][:cyclic] = @constraint(mod, [jd=JD, jd_a=JD_A],
        SOC_AD_0[first(JD_A)] <= SOC_AD_0[last(JD_A)] + sum(V[last(JD_A),jd]*delta_e[jd] for jd in JD)  )

        # Inspired by me 
        #### check all the initial opt problem constraint (in particular 21,22,23,24), maybe define SOC_repr_0

        # Charge and discharge limits - (19),(20)
        mod.ext[:constraints][:dischargin_lim] = @constraint(mod, [jh=JH,jd=JD], 
        dhH[jh,jd] <= capH/1000 )  # [TWh]
        
        mod.ext[:constraints][:charging_lim] = @constraint(mod, [jh=JH,jd=JD], 
        chH[jh,jd] <= capH/1000 ) # [TWh]

        # SOC upper limit - (23),(24) - maybe redundant, but before is done only for SOC_AD
        mod.ext[:constraints][:SOC_upper_limit_repr] = @constraint(mod, [jh = JH, jd = JD],
        SOC[jh,jd] <= volH)

        mod.ext[:constraints][:SOC_upper_limit_repr_0] = @constraint(mod, [jd=JD],
        SOC_0[jd] <= volH)

        # cyclic balance - divided into two constraints to consider the first hour of the day that has to refer to the last one of the previous day 
        # @constraint(mod, [jh=first(JH),jd=JD[2:1:end]], SOC[jh,jd] == SOC[last(JH),jd-1] + η_ch*chH[jh,jd] - dhH[jh,jd]/η_dh ) # short term storage

        ## SOC dynamics and balance
        # SOC dynamics for representative periods (intraperiod) - (21),(22)


        mod.ext[:constraints][:SOC_update_intra_0] = @constraint(mod, [jh=first(JH),jd=JD],  
        SOC[jh,jd] == SOC_0[jd] + sum(η_ch*chH[jh,jd] - dhH[jh,jd]/η_dh) )

        mod.ext[:constraints][:SOC_update_intra] = @constraint(mod, [jh=JH[2:1:end],jd=JD], 
        SOC[jh, jd] == SOC[jh-1, jd] + (η_ch * chH[jh, jd] - dhH[jh, jd] / η_dh) )         # *1h is implied

        # SOC dynamics interperiod (maybe useful for plots later) - TO BE CHANGED EVENTUALLY
        # mod.ext[:constraints][:SOC_update_inter0] = @constraint(mod, [jh=first(JH), jd_a=JD_A],
        # SOC[jh,jd_a] = SOC[last(JH),jd_a-1] + sum(chron[jd,jd_a]* (η_ch*ch[jh,jd] - dhH[jh,jd]/η_dh)
        # for jd in JD)) 
        # mod.ext[:constraints][:SOC_update_inter] = @constraint(mod, [jh=JH[2:end], jd_a=JD_A],
        # SOC[jh,jd_a] = SOC[jh-1,jd_a] + sum(chron[jd,jd_a]* (η_ch*ch[jh,jd] - dhH[jh,jd]/η_dh)
        # for jd in JD))
        # the last two can be use with P2AP funtion
    end
    
    return mod

end