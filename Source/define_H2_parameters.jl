function define_H2_parameters!(H2::Dict,data::Dict,ts::DataFrame,repr_days::DataFrame)

    # H2 demand
    H2["D"] = [data["conv_factor"]*data["demand"]/8760 for h in 1:data["nTimesteps"], d in 1:data["nReprDays"]] 

    # Weights
    H2["W"] = W = Dict(jd => repr_days[!,:weights][jd] for jd=1:data["nReprDays"])

    return H2
end