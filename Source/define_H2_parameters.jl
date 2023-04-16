function define_H2_parameters!(H2::Dict, data::Dict, ts::Dict, repr_days::Dict)
    
    H2["D"] = Array{Float64}(undef, data["nTimesteps"], data["nReprDays"],data["nYears"])
    H2["elasticity"] = Array{Float64}(undef, data["nYears"])
    H2["W"] = Dict()
    
    for jy in keys(years)
        # Timeseries
        # H2["D"] = [data["conv_factor"]*data["demand"]/8760 for h in 1:data["nTimesteps"], d in 1:data["nReprDays"]] # TWh/h
        H2["D"][:,:,jy] = [ts[jy][!,:LOAD_H2][round(Int, data["nTimesteps"]*(repr_days[jy][!,:periods][jd]-1) + jh)]/10^6 for jh=1:data["nTimesteps"], jd=1:data["nReprDays"]] # from MWh to TWh
        # N.B. to change H2 load in "timeseries.csv", go back to "H2demandincreas" excel sheet

        H2["elasticity"][jy] = data["elastic_sec_H2"]

        # Weights
        H2["W"][jy] = W = Dict(jd => repr_days[jy][!,:weights][jd] for jd=1:data["nReprDays"])
    end

    return H2
end