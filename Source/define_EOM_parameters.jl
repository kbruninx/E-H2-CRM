function define_EOM_parameters!(EOM::Dict, data::Dict, ts::Dict, repr_days::Dict)
    
    EOM["D"] = Array{Float64}(undef, data["nTimesteps"], data["nReprDays"],data["nYears"])          # energy demand
    EOM["elasticity"] = Array{Float64}(undef, data["nTimesteps"], data["nReprDays"],data["nYears"])
    EOM["W"] = Dict()

    for jy in keys(years)
        # timeseries
        EOM["D"][:,:,jy] = [ts[jy][!,:LOAD][round(Int, data["nTimesteps"]*(repr_days[jy][!,:periods][jd]-1) + jh)]/10^5 for jh=1:data["nTimesteps"], jd=1:data["nReprDays"]] # 10^2 GWh
    
        EOM["elasticity"][:,:,jy] = [ts[jy][!,:ELASTICITY_EL][round(Int, data["nTimesteps"]*(repr_days[jy][!,:periods][jd]-1) + jh)]/10^5 for jh=1:data["nTimesteps"], jd=1:data["nReprDays"]] # 10^2 GWh

        # weights of representative days
        EOM["W"][jy] = W = Dict(jd => repr_days[jy][!,:weights][jd] for jd=1:data["nReprDays"])
    end

    return EOM
end