function build_h2s_agent!(m::String, mod::Model, H2::Dict)
    # Extract sets
    JH = mod.ext[:sets][:JH]
    JD = mod.ext[:sets][:JD]
    JY = mod.ext[:sets][:JY]

    # Extract common parameters
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
    gH_bar = mod.ext[:parameters][:gH_bar] # element in ADMM penalty term related to hydrogen market
    ρ_H2 = mod.ext[:parameters][:ρ_H2] # rho-value in ADMM related to H2 market
    λ_CM = mod.ext[:parameters][:λ_CM] # capacity market prices
    cap_bar = mod.ext[:parameters][:cap_bar]  # element in ADMM penalty term related to CM
    ρ_CM = mod.ext[:parameters][:ρ_CM]  # rho-value in ADMM related to CM auctions
    λ_HCM = mod.ext[:parameters][:λ_HCM] # hydrogen capacity market prices
    capH_bar = mod.ext[:parameters][:capH_bar]  # element in ADMM penalty term related to hydrogen CM
    ρ_HCM = mod.ext[:parameters][:ρ_HCM]  # rho-value in ADMM related to hydrogen CM auctions


    if m == "H2demand"
        # Extract parameters
        WTP = mod.ext[:parameters][:WTP]  # willingness to pay ("price cap") of consumers
        ela_H2 = H2["elasticity"] # section of price-elasticity [MWh]

        # Create variables
        gH_VOLL = mod.ext[:variables][:gH_VOLL] = @variable(mod, [jh = JH, jd = JD, jy = JY], lower_bound = 0, base_name = "H2_inelastic_demand")
        gH_ela = mod.ext[:variables][:gH_ela] = @variable(mod, [jh = JH, jd = JD, jy = JY], lower_bound = 0, base_name = "H2_elastic_demand")
        α = mod.ext[:variables][:α] = @variable(mod, base_name = "VAR")
        u = mod.ext[:variables][:u] = @variable(mod, [jy = JY], lower_bound = 0, base_name = "tail profit difference")  # profit difference of worst-case tail scenarios with respect to the VAR

        # Warm start
        #=
        for jy in JY
            start_value = Matrix(CSV.read(joinpath(home_dir, "Input", "RA10_elCM", string("Elastic_demand_H2_", jy, ".csv")), delim=";", DataFrame))
            for jh in JH, jd in JD, jy in JY
                set_start_value(gH_ela[jh, jd, jy], start_value[jh, jd])
            end
        end
        =#

        # Create expressions

        gH_positive = mod.ext[:expressions][:gH_positive] = @expression(mod, gH_VOLL + gH_ela)
        gH = mod.ext[:expressions][:gH] = @expression(mod, - gH_positive)

        profit = mod.ext[:expressions][:profit] = @expression(mod, [jy = JY],
        sum(W[jd, jy] * (WTP * gH_positive[jh, jd, jy] - (gH_ela[jh, jd, jy])^2 * WTP / (2 * ela_H2[jy])
        - λ_H2[jh, jd, jy] * gH_positive[jh, jd, jy]) for jh in JH, jd in JD)
        - σH * λ_HCM * HCM["D"])

        cost = mod.ext[:expressions][:cost] = @expression(mod, [jy = JY],
        sum(W[jd, jy] * λ_H2[jh, jd, jy] * gH_positive[jh, jd, jy] for jh in JH, jd in JD)
        + σH * λ_HCM * HCM["D"])

        CVAR = mod.ext[:expressions][:CVAR] = @expression(mod,
        α - ((1 / β) * sum(P[jy] * u[jy] for jy in JY)))

        # Objective 
        mod.ext[:objective] = @objective(mod, Min,
        - γ * sum(P[jy] * profit[jy] for jy in JY) - (1 - γ) * CVAR
        + sum(ρ_H2 / 2 * W[jd, jy] * (gH[jh, jd, jy] - gH_bar[jh, jd, jy])^2 for jh in JH, jd in JD, jy in JY)
        ) 

        # Constraints
        mod.ext[:constraints][:elastic_demand_limit] = @constraint(mod, [jh = JH, jd = JD, jy = JY],
        gH_ela[jh, jd, jy] <= ela_H2[jy])

        mod.ext[:constraints][:demand_limit] = @constraint(mod, [jh = JH, jd = JD, jy = JY],
        gH_VOLL[jh, jd, jy] <= 0.8 * H2["D"][jh, jd, jy])

        if γ < 1
            # CVAR constraint
            mod.ext[:constraints][:VAR_threshold] = @constraint(mod, [jy = JY],
            α - profit[jy] <= u[jy])        
        end

    elseif m == "electrolysis"

        # Extract parameters
        η_E_H2 = mod.ext[:parameters][:η_E_H2] # efficiency electrolysis
        IC = mod.ext[:parameters][:IC] # annuity investment costs

        # Decision variables
        capH = mod.ext[:variables][:capH] = @variable(mod, lower_bound = 0, base_name = "capacity")
        capH_cm = mod.ext[:variables][:capH_cm] = @variable(mod, lower_bound = 0, base_name = "capacity offered HCM")
        g = mod.ext[:variables][:g] = @variable(mod, [jh = JH, jd = JD, jy = JY], upper_bound = 0, base_name = "demand_electricity_hydrogen") # note this is defined as a negative number, consumption
        α = mod.ext[:variables][:α] = @variable(mod, base_name = "VAR")
        u = mod.ext[:variables][:u] = @variable(mod, [jy = JY], lower_bound = 0, base_name = "tail profit difference")  # profit difference of worst-case tail scenarios with respect to the VAR

        # Warm start
        #=
        for jy in JY
            start_value = Matrix(CSV.read(joinpath(home_dir, "Input", "RA10_elCM", string("H2Generationelectrolysis_", jy, ".csv")), delim=";", DataFrame)) ./ (-η_E_H2)
            for jh in JH, jd in JD, jy in JY
                set_start_value(g[jh, jd, jy], start_value[jh, jd])
            end
        end
        set_start_value(capH, (CSV.read(joinpath(home_dir, "Input", "RA10_elCM", "capacity_h2s.csv"), delim=";", DataFrame))[!, :CAP_electrolysis][1])
        =#

        # Create affine expressions  
        #inv_cost = mod.ext[:expressions][:inv_cost] = @expression(mod, IC * capH)
        
        gH = mod.ext[:expressions][:gH] = @expression(mod, -η_E_H2 * g) # [TWh]

        cost = mod.ext[:expressions][:cost] = @expression(mod, [jy = JY],
        IC * capH
        - sum(W[jd, jy] * λ_EOM[jh, jd, jy] * g[jh, jd, jy] for jh in JH, jd in JD))

        revenue = mod.ext[:expressions][:revenue] = @expression(mod, [jy = JY],
        sum(W[jd, jy] * λ_H2[jh, jd, jy] * gH[jh, jd, jy] for jh in JH, jd in JD)
        + σH * λ_HCM * capH_cm)

        profit = mod.ext[:expressions][:profit] = @expression(mod, [jy = JY],
        revenue[jy] - cost[jy] )

        CVAR = mod.ext[:expressions][:CVAR] = @expression(mod,
        α - ((1 / β) * sum(P[jy] * u[jy] for jy in JY)))

        # Definition of the objective function
        mod.ext[:objective] = @objective(mod, Min,
        - γ * sum(P[jy] * (revenue[jy] - cost[jy]) for jy in JY) - (1 - γ) * CVAR
        + sum(ρ_EOM / 2 * W[jd, jy] * (g[jh, jd, jy] - g_bar[jh, jd, jy])^2 for jh in JH, jd in JD, jy in JY)
        + sum(ρ_H2 / 2 * W[jd, jy] * (gH[jh, jd, jy] - gH_bar[jh, jd, jy])^2 for jh in JH, jd in JD, jy in JY)
        + σH * ρ_HCM / 2 * (capH_cm - capH_bar)^2
        )

        # Constraints

        # Electricity consumption
        mod.ext[:constraints][:elec_consumption] = @constraint(mod, [jh = JH, jd = JD, jy = JY],
        -g[jh, jd, jy] <= capH /100  # [10^2 GWh]
        )

        # Capacity offered in H2 CM - they offer up to their installed capacity (corrected by thir efficiency to represent true output)
        if σH == 1
            mod.ext[:constraints][:cap_offer] = @constraint(mod,
            capH_cm <= capH * η_E_H2)  # previously no * η_E_H2, since the capacity demand refers alreay to these
        end

        if γ < 1
            # CVAR constraint
            mod.ext[:constraints][:VAR_threshold] = @constraint(mod, [jy = JY],
            α - (revenue[jy] - cost[jy]) <= u[jy])
        end

    elseif m == "H2storage"

        # Extract parameters
        JD_A = mod.ext[:sets][:JD_A] # sets of all days used for storage SOC
        index_repr = mod.ext[:sets][:index_repr] # index of representative days
        V = mod.ext[:sets][:order_matrix]
        VC = mod.ext[:parameters][:VC] # variable cost storage
        η_ch = mod.ext[:parameters][:η_ch] # charging efficiency 
        η_dh = mod.ext[:parameters][:η_dh] # discharging efficiency
        IC_cap = mod.ext[:parameters][:IC_cap] # annuity investment cost for capacity
        IC_vol = mod.ext[:parameters][:IC_vol] # annuity investment cost for volume

        # Decision variables
        capH = mod.ext[:variables][:capH] = @variable(mod, lower_bound = 0, base_name = "capacity") # MW
        volH = mod.ext[:variables][:volH] = @variable(mod, lower_bound = 0, base_name = "volume") # MWh
        capH_cm = mod.ext[:variables][:capH_cm] = @variable(mod, lower_bound = 0, base_name = "capacity offered HCM") # MW

        chH = mod.ext[:variables][:chH] = @variable(mod, [jh = JH, jd = JD, jy = JY], lower_bound = 0, base_name = "hydrogen_charging")
        dhH = mod.ext[:variables][:dhH] = @variable(mod, [jh = JH, jd = JD, jy = JY], lower_bound = 0, base_name = "hydrogen_discharging")
        delta_e = mod.ext[:expressions][:delta_e] = @variable(mod, [jd = JD, jy = JY])
        SOC = mod.ext[:variables][:SOC] = @variable(mod, [jh = JH, jd = JD, jy = JY], lower_bound = 0, base_name = "state_of_charge")
        SOC_0 = mod.ext[:variables][:SOC_0] = @variable(mod, [jd = JD, jy = JY], lower_bound = 0, base_name = "initial_state_of_charge_repr_day")    # state of charge at the start of a period
        SOC_AD_0 = mod.ext[:variables][:SOC_AD_0] = @variable(mod, [jd_a = JD_A, jy = JY], lower_bound = 0, base_name = "initial_state_of_charge_all_days")
        max_deviation_SOC = mod.ext[:variables][:max_deviation_SOC] = @variable(mod, [jd = JD, jy = JY], lower_bound = 0, base_name = "max_deviaion_SOC")
        min_deviation_SOC = mod.ext[:variables][:min_deviation_SOC] = @variable(mod, [jd = JD, jy = JY], lower_bound = 0, base_name = "min_deviaion_SOC")    # refer to a negative deviation
        max_deviation_SOC_AD = mod.ext[:variables][:max_deviation_SOC_AD] = @variable(mod, [jd_a = JD_A, jy = JY], lower_bound = 0)
        min_deviation_SOC_AD = mod.ext[:variables][:min_deviation_SOC_AD] = @variable(mod, [jd_a = JD_A, jy = JY], lower_bound = 0)
        α = mod.ext[:variables][:α] = @variable(mod, base_name = "VAR")
        u = mod.ext[:variables][:u] = @variable(mod, [jy = JY], lower_bound = 0, base_name = "tail profit difference")  # profit difference of worst-case tail scenarios with respect to the VAR

        # Warm start
        #=
        for jy in JY
            start_value_chH = Matrix(CSV.read(joinpath(home_dir, "Input", "RA10_elCM", string("ChargeH2storage_", jy, ".csv")), delim=";", DataFrame))
            start_value_dhH = Matrix(CSV.read(joinpath(home_dir, "Input", "RA10_elCM", string("DischargeH2storage_", jy, ".csv")), delim=";", DataFrame))
            for jh in JH, jd in JD, jy in JY
                set_start_value(chH[jh, jd, jy], start_value_chH[jh, jd])
                set_start_value(dhH[jh, jd, jy], start_value_dhH[jh, jd])
            end
        end
        set_start_value(capH, (CSV.read(joinpath(home_dir, "Input", "RA10_elCM", "capacity_h2s.csv"), delim=";", DataFrame))[!, :CAP_H2storage][1])
        set_start_value(volH, (CSV.read(joinpath(home_dir, "Input", "RA10_elCM", "storage_volume.csv"), delim=";", DataFrame))[!, :x1][1])
        =#
        
        # Create affine expressions  
        #inv_cost = mod.ext[:expressions][:inv_cost] = @expression(mod, IC_cap * capH + IC_vol * volH)
        gH = mod.ext[:expressions][:gH] = @expression(mod, dhH - chH)

        cost = mod.ext[:expressions][:cost] = @expression(mod, [jy = JY],
        IC_cap * capH + IC_vol * volH
        + sum(W[jd, jy] * λ_H2[jh, jd, jy] * chH[jh, jd, jy] for jh in JH, jd in JD))

        revenue = mod.ext[:expressions][:revenue] = @expression(mod, [jy = JY],
        sum(W[jd, jy] * λ_H2[jh, jd, jy] * dhH[jh, jd, jy] for jh in JH, jd in JD)
        + σH * λ_HCM * capH_cm
        )

        profit = mod.ext[:expressions][:profit] = @expression(mod, [jy = JY],
        revenue[jy] - cost[jy] )

        CVAR = mod.ext[:expressions][:CVAR] = @expression(mod,
        α - ((1 / β) * sum(P[jy] * u[jy] for jy in JY)))

        # Definition of the objective function
        mod.ext[:objective] = @objective(mod, Min,
        - γ * sum(P[jy] * (revenue[jy] - cost[jy]) for jy in JY) - (1 - γ) * CVAR
        + sum(ρ_H2 / 2 * W[jd, jy] * (gH[jh, jd, jy] - gH_bar[jh, jd, jy])^2 for jh in JH, jd in JD, jy in JY)
        + σH * ρ_HCM / 2 * (capH_cm - capH_bar)^2
        )

        # Constraints

        # Inspired by Gonzato, Bruninx and Delarue - Long term storage in generation expansion planning models with a reduced temporal scope
        #### take in input V (ordering_variable), check dimensions of V in the constraints

        # Linking SOC_0 and SOC_AD_0 
        mod.ext[:constraints][:D2AD] = @constraint(mod, [jd = JD, jd_a = JD_A, jy = JY],
            SOC_0[jd, jy] == SOC_AD_0[Int(index_repr[jd]), jy])

        # Define energy change for each representative day - (33)
        mod.ext[:constraints][:delta_e_repr] = @constraint(mod, [jh = JH, jd = JD, jy = JY],
            delta_e[jd, jy] == sum((η_ch * chH[jh, jd, jy] - dhH[jh, jd, jy] / η_dh) for jh in JH))

        # Base state of charge definition (for all days) - (34) - shouldn't they be of the previous jd_a ?
        mod.ext[:constraints][:SOC_AD_0] = @constraint(mod, [jd = JD, jd_a = JD_A[2:end], jy = JY],
            SOC_AD_0[jd_a, jy] == SOC_AD_0[jd_a-1, jy] + sum(V[jd_a, jd, jy] * delta_e[jd, jy] for jd in JD))

        # Max positive and negative deviations from the base state of charge for representiative periods  - (35),(36)
        mod.ext[:constraints][:max_deviation_repr] = @constraint(mod, [jh = JH, jd = JD, jy = JY],
            max_deviation_SOC[jd, jy] >= SOC[jh, jd, jy] - SOC_0[jd, jy])

        mod.ext[:constraints][:min_deviation_repr] = @constraint(mod, [jh = JH, jd = JD, jy = JY],
            min_deviation_SOC[jd, jy] >= SOC_0[jd, jy] - SOC[jh, jd, jy])

        # Max positive and negative deviations from the base state of charge for non-representiative periods - (37),(38)
        mod.ext[:constraints][:max_deviation_all] = @constraint(mod, [jd = JD, jd_a = JD_A, jy = JY],
            max_deviation_SOC_AD[jd_a, jy] == sum(V[jd_a, jd, jy] * max_deviation_SOC[jd, jy] for jd in JD))

        mod.ext[:constraints][:min_deviation_all] = @constraint(mod, [jd = JD, jd_a = JD_A, jy = JY],
            min_deviation_SOC_AD[jd_a, jy] == sum(V[jd_a, jd, jy] * min_deviation_SOC[jd, jy] for jd in JD))

        # SOC limits for all periods - (39),(40)
        mod.ext[:constraints][:SOC_upper_limit_all] = @constraint(mod, [jd_a = JD_A, jy = JY],
            volH >= SOC_AD_0[jd_a, jy] + max_deviation_SOC_AD[jd_a, jy])

        mod.ext[:constraints][:SOC_lower_limit_all] = @constraint(mod, [jd_a = JD_A, jy = JY],
            SOC_AD_0[jd_a, jy] - min_deviation_SOC_AD[jd_a, jy] >= 0)

        # Cyclic constraint - (41) but adpated, the one in the paper does not make sense
        mod.ext[:constraints][:cyclic] = @constraint(mod, [jd = JD, jd_a = JD_A, jy = JY],
            SOC_AD_0[first(JD_A), jy] <= SOC_AD_0[last(JD_A), jy] + sum(V[last(JD_A), jd, jy] * delta_e[jd, jy] for jd in JD))

        # Inspired by me 
        #### check all the initial opt problem constraint (in particular 21,22,23,24), maybe define SOC_repr_0

        # Charge and discharge limits - (19),(20)
        mod.ext[:constraints][:dischargin_lim] = @constraint(mod, [jh = JH, jd = JD, jy = JY],
            dhH[jh, jd, jy] <= capH / 100)  # [10^2 GWh]

        mod.ext[:constraints][:charging_lim] = @constraint(mod, [jh = JH, jd = JD, jy = JY],
            chH[jh, jd, jy] <= capH / 100) # [10^2 GWh]

        # SOC upper limit - (23),(24) - maybe redundant, but before is done only for SOC_AD
        mod.ext[:constraints][:SOC_upper_limit_repr] = @constraint(mod, [jh = JH, jd = JD, jy = JY],
            SOC[jh, jd, jy] <= volH)

        mod.ext[:constraints][:SOC_upper_limit_repr_0] = @constraint(mod, [jd = JD, jy = JY],
            SOC_0[jd, jy] <= volH)

        # cyclic balance - divided into two constraints to consider the first hour of the day that has to refer to the last one of the previous day 
        # @constraint(mod, [jh=first(JH),jd=JD[2:1:end]], SOC[jh,jd] == SOC[last(JH),jd-1] + η_ch*chH[jh,jd] - dhH[jh,jd]/η_dh ) # short term storage

        ## SOC dynamics and balance
        # SOC dynamics for representative periods (intraperiod) - (21),(22)


        mod.ext[:constraints][:SOC_update_intra_0] = @constraint(mod, [jh = first(JH), jd = JD, jy = JY],
            SOC[jh, jd, jy] == SOC_0[jd, jy] + (η_ch * chH[jh, jd, jy] - dhH[jh, jd, jy] / η_dh))

        mod.ext[:constraints][:SOC_update_intra] = @constraint(mod, [jh = JH[2:1:end], jd = JD, jy = JY],
            SOC[jh, jd, jy] == SOC[jh-1, jd, jy] + (η_ch * chH[jh, jd, jy] - dhH[jh, jd, jy] / η_dh))         # *1h is implied


        # SOC dynamics interperiod (maybe useful for plots later) - TO BE CHANGED EVENTUALLY
        # mod.ext[:constraints][:SOC_update_inter0] = @constraint(mod, [jh=first(JH), jd_a=JD_A],
        # SOC[jh,jd_a] = SOC[last(JH),jd_a-1] + sum(chron[jd,jd_a]* (η_ch*ch[jh,jd] - dhH[jh,jd]/η_dh)
        # for jd in JD)) 
        # mod.ext[:constraints][:SOC_update_inter] = @constraint(mod, [jh=JH[2:end], jd_a=JD_A],
        # SOC[jh,jd_a] = SOC[jh-1,jd_a] + sum(chron[jd,jd_a]* (η_ch*ch[jh,jd] - dhH[jh,jd]/η_dh)
        # for jd in JD))
        # the last two can be use with P2AP funtion

        if γ < 1
            # CVAR constraint
            mod.ext[:constraints][:VAR_threshold] = @constraint(mod, [jy = JY],
            α - (revenue[jy] - cost[jy]) <= u[jy])
        end

        if σH == 1   # capacity offered in the H2 CM (correctd for the discharge efficiency)
            mod.ext[:constraints][:cap_offer] = @constraint(mod,
            capH_cm <= capH)
        end
    end

    return mod

end