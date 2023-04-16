function define_ps_parameters!(m::String, mod::Model, data::Dict, ts::Dict, repr_days::Dict)
    # Parameters 
    if m == "Edemand"
        mod.ext[:parameters][:WTP] = data["WTP"]
    else
        mod.ext[:parameters][:VC] = data["VC"]
        mod.ext[:parameters][:IC] = data["OC"]/data["Lifetime"] * (1 + data["inflation"])^data["Lifetime"] / (1 + data["discount_rate"])^data["Lifetime"]

        if m == "H2turbine"
            mod.ext[:parameters][:η_H2_E] = data["efficiency_H2_E"]
        end

        # Availability factors
        if data["AF_ts"] != "NA" 
            mod.ext[:timeseries][:AF] = zeros(data["nTimesteps"],data["nReprDays"],data["nYears"])
            mod.ext[:timeseries][:AF] = [ts[jy][!,data["AF_ts"]][round(Int,data["nTimesteps"]*repr_days[jy][!,:periods][jd]+jh)] for jh=1:data["nTimesteps"], jd=1:data["nReprDays"], jy=1:data["nYears"]]  
        else
            mod.ext[:timeseries][:AF] = ones(data["nTimesteps"],data["nReprDays"],data["nYears"]) 
        end
    end
   
    return mod
end