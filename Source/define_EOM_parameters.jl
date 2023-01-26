function define_EOM_parameters!(EOM::Dict,data::Dict,ts::DataFrame,repr_days::DataFrame)
    # timeseries
    EOM["D"] = [ts[!,:LOAD][round(Int,data["nTimesteps"]*repr_days[!,:periods][jd]+jh)]/10^6 for jh=1:data["nTimesteps"], jd=1:data["nReprDays"]] # from MWh to TWh

    # weights of representative days
    EOM["W"] = W = Dict(jd => repr_days[!,:weights][jd] for jd=1:data["nReprDays"])

    return EOM
end