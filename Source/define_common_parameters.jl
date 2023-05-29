function define_common_parameters!(m::String, mod::Model, data::Dict, ts::Dict, repr_days::Dict, agents::Dict)
    # Solver settings
    # Define dictonaries for sets, parameters, timeseries, variables, constraints & expressions
    mod.ext[:sets] = Dict()
    mod.ext[:parameters] = Dict()
    mod.ext[:timeseries] = Dict()
    mod.ext[:variables] = Dict()
    mod.ext[:constraints] = Dict()
    mod.ext[:expressions] = Dict()

    # Sets
    mod.ext[:sets][:JY] = 1:data["nYears"]
    mod.ext[:sets][:JD] = 1:data["nReprDays"]
    mod.ext[:sets][:JH] = 1:data["nTimesteps"]

    # Parameters
    mod.ext[:parameters][:W] = [repr_days[jy][!, :weights][jd] for jd in mod.ext[:sets][:JD], jy in mod.ext[:sets][:JY]]      # weights of each representative day
    mod.ext[:parameters][:P] = ones(data["nYears"]) / data["nYears"]    # probability of each scenario - uniform distribution
    mod.ext[:parameters][:γ] = data["gamma"]            # weight of expected revenues and CVAR
    mod.ext[:parameters][:β] = data["beta"]          # risk aversion parameter - represents the cumulative probability of worst-case scenarios
    mod.ext[:parameters][:σ] = data["sigmaCM"] # switch for the capacity market
    mod.ext[:parameters][:σH] = data["sigmaHCM"] # switch for the capacity market

    # Parameters related to the EOM
    mod.ext[:parameters][:λ_EOM] = zeros(data["nTimesteps"], data["nReprDays"], data["nYears"])   # Price structure
    mod.ext[:parameters][:g_bar] = zeros(data["nTimesteps"], data["nReprDays"], data["nYears"])   # ADMM penalty term
    mod.ext[:parameters][:ρ_EOM] = data["rho_EOM"]                                                   # ADMM rho value 

    # Parameters related to the H2 market
    mod.ext[:parameters][:λ_H2] = zeros(data["nTimesteps"], data["nReprDays"], data["nYears"])       # Price structure
    mod.ext[:parameters][:gH_bar] = zeros(data["nTimesteps"], data["nReprDays"], data["nYears"])     # ADMM penalty term
    mod.ext[:parameters][:ρ_H2] = data["rho_H2"]                                           # ADMM rho value 

    # Parameters related to the electricity CM
    mod.ext[:parameters][:λ_CM] = 0    # Price structure
    mod.ext[:parameters][:cap_bar] = 0   # ADMM penalty term
    mod.ext[:parameters][:ρ_CM] = data["rho_CM"]           # ADMM rho value 

    # Parameters related to the hydrogen CM
    mod.ext[:parameters][:λ_HCM] = 0    # Price structure
    mod.ext[:parameters][:capH_bar] = 0   # ADMM penalty term
    mod.ext[:parameters][:ρ_HCM] = data["rho_HCM"]           # ADMM rho value 

    # Covered by EOM?
    if data["EOM"] == "YES"
        mod.ext[:parameters][:EOM] = 1
        push!(agents[:eom], m)
    else
        mod.ext[:parameters][:EOM] = 0
    end

    # Covered by Hydrogen Market 
    if data["H2"] == "YES"
        mod.ext[:parameters][:H2] = 1
        push!(agents[:h2], m)
    else
        mod.ext[:parameters][:H2] = 0
    end

    # Covered by Electrcity Capacity Market 
    if data["CM"] == "YES"
        mod.ext[:parameters][:CM] = 1
        push!(agents[:cm], m)
    else
        mod.ext[:parameters][:CM] = 0
    end

    # Covered by Hydrogen Capacity Market 
    if data["HCM"] == "YES"
        mod.ext[:parameters][:HCM] = 1
        push!(agents[:hcm], m)
    else
        mod.ext[:parameters][:HCM] = 0
    end

    return mod, agents
end