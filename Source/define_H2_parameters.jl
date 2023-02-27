function define_H2_parameters!(H2::Dict,data::Dict,ts::DataFrame,repr_days::DataFrame)

    # H2 demand
    # H2["D"] = [data["conv_factor"]*data["demand"]/8760 for h in 1:data["nTimesteps"], d in 1:data["nReprDays"]] # TWh/h
    H2["D"] = [ts[!,:LOAD_H2][round(Int,data["nTimesteps"]*repr_days[!,:periods][jd]+jh)]/10^6 for jh=1:data["nTimesteps"], jd=1:data["nReprDays"]] # from MWh to TWh
    # N.B. to change H2 load in "timeseries.csv", go back to "H2demandincreas" excel sheet

    H2["elasticity"] = data["elasticity_H2"]

    # Weights
    H2["W"] = W = Dict(jd => repr_days[!,:weights][jd] for jd=1:data["nReprDays"])

    return H2
end