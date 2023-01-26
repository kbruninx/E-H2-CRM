function define_common_parameters!(m::String,mod::Model, data::Dict, ts::DataFrame, repr_days::DataFrame, agents::Dict)
    # Solver settings
    # Define dictonaries for sets, parameters, timeseries, variables, constraints & expressions
    mod.ext[:sets] = Dict()
    mod.ext[:parameters] = Dict()
    mod.ext[:timeseries] = Dict()
    mod.ext[:variables] = Dict()
    mod.ext[:constraints] = Dict()
    mod.ext[:expressions] = Dict()

    # Sets
    mod.ext[:sets][:JD] = 1:data["nReprDays"]
    mod.ext[:sets][:JH] = 1:data["nTimesteps"]

    # Parameters
    mod.ext[:parameters][:W] = [repr_days[!,:weights][jd] for jd in mod.ext[:sets][:JD]]      # weights of each representative day

    # Parameters related to the EOM
    mod.ext[:parameters][:λ_EOM] = zeros(data["nTimesteps"],data["nReprDays"])   # Price structure
    mod.ext[:parameters][:g_bar] = zeros(data["nTimesteps"],data["nReprDays"])   # ADMM penalty term
    mod.ext[:parameters][:ρ_EOM] = data["rho_EOM"]                               # ADMM rho value 
    
    # Parameters related to the H2 market
    mod.ext[:parameters][:λ_H2] = zeros(data["nTimesteps"],data["nReprDays"])       # Price structure
    mod.ext[:parameters][:gH_bar] = zeros(data["nTimesteps"],data["nReprDays"])     # ADMM penalty term
    mod.ext[:parameters][:ρ_H2] = data["rho_H2"]                                    # ADMM rho value 
    
    # Covered by EOM?
    if data["EOM"] == "YES" 
        mod.ext[:parameters][:EOM] = 1
        push!(agents[:eom],m)
    else
        mod.ext[:parameters][:EOM] = 0
    end

    # Covered by Hydrogen Market 
    if data["H2"] == "YES" 
        mod.ext[:parameters][:H2] = 1
        push!(agents[:h2],m)
    else
        mod.ext[:parameters][:H2] = 0
    end

    return mod, agents
end