function define_ps_parameters!(mod::Model, data::Dict,ts::DataFrame,repr_days::DataFrame)
    # Parameters 
    mod.ext[:parameters][:VC] = data["VC"]
    mod.ext[:parameters][:IC] = data["OC"]/data["Lifetime"]
    if mod.ext[:parameters][:EOM] == 1 && mod.ext[:parameters][:H2] == 1
        mod.ext[:parameters][:η_H2_E] = data["efficiency_H2_E"]
    end
    
   # Availability factors
    if data["AF_ts"] != "NA" 
        mod.ext[:timeseries][:AF] = [ts[!,data["AF_ts"]][round(Int,data["nTimesteps"]*repr_days[!,:periods][jd]+jh)] for jh=1:data["nTimesteps"], jd=1:data["nReprDays"]]  
    else
        mod.ext[:timeseries][:AF] = ones(data["nTimesteps"],data["nReprDays"]) 
    end 
   
    return mod
end