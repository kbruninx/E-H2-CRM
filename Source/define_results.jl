function define_results!(data::Dict, results::Dict, ADMM::Dict, agents::Dict, EOM::Dict, H2::Dict)
    results["g"] = Dict()
    results["h2"] = Dict()
    results["cap_cm"] = Dict()
    results["capH_cm"] = Dict()
    results["SOC"] = Dict() # probably not needed

    for m in agents[:eom]
        results["g"][m] = CircularBuffer{Array{Float64,3}}(data["CircularBufferSize"])
        push!(results["g"][m], zeros(data["nTimesteps"], data["nReprDays"], data["nYears"]))
        if m == "Edemand"
            results["g_ela"] = CircularBuffer{Array{Float64,3}}(data["CircularBufferSize"])
            push!(results["g_ela"], zeros(data["nTimesteps"], data["nReprDays"], data["nYears"]))
        end
    end
    for m in agents[:h2]
        results["h2"][m] = CircularBuffer{Array{Float64,3}}(data["CircularBufferSize"])
        push!(results["h2"][m], zeros(data["nTimesteps"], data["nReprDays"], data["nYears"]))
        if m == "H2storage"
            results["SOC"] = CircularBuffer{Array{Float64,3}}(data["CircularBufferSize"])
            push!(results["SOC"], zeros(data["nTimesteps"], data["nReprDays"], data["nYears"]))
            results["SOC_AD_0"] = CircularBuffer{Array{Float64,2}}(data["CircularBufferSize"])
            push!(results["SOC_AD_0"], zeros(data["nDays"], data["nYears"]))
            results["dhH"] = CircularBuffer{Array{Float64,3}}(data["CircularBufferSize"])
            push!(results["dhH"], zeros(data["nTimesteps"], data["nReprDays"], data["nYears"]))
            results["chH"] = CircularBuffer{Array{Float64,3}}(data["CircularBufferSize"])
            push!(results["chH"], zeros(data["nTimesteps"], data["nReprDays"], data["nYears"]))
        end
        if m == "H2demand"
            results["gH2_ela"] = CircularBuffer{Array{Float64,3}}(data["CircularBufferSize"])
            push!(results["gH2_ela"], zeros(data["nTimesteps"], data["nReprDays"], data["nYears"]))
        end
    end
    for m in agents[:cm]
        results["cap_cm"][m] = CircularBuffer{Float64}(data["CircularBufferSize"])
        push!(results["cap_cm"][m], 0.0)
    end
    for m in agents[:hcm]
        results["capH_cm"][m] = CircularBuffer{Float64}(data["CircularBufferSize"])
        push!(results["capH_cm"][m], 0.0)
    end



    results["λ"] = Dict()
    results["λ"]["EOM"] = CircularBuffer{Array{Float64,3}}(data["CircularBufferSize"])
    push!(results["λ"]["EOM"], zeros(data["nTimesteps"], data["nReprDays"], data["nYears"]))
    results["λ"]["H2"] = CircularBuffer{Array{Float64,3}}(data["CircularBufferSize"])
    push!(results["λ"]["H2"], zeros(data["nTimesteps"], data["nReprDays"], data["nYears"]))
    results["λ"]["CM"] = CircularBuffer{Float64}(data["CircularBufferSize"])
    push!(results["λ"]["CM"], 0.0)
    results["λ"]["HCM"] = CircularBuffer{Float64}(data["CircularBufferSize"])
    push!(results["λ"]["HCM"], 0.0)

    ADMM["Imbalances"] = Dict()
    ADMM["Imbalances"]["EOM"] = CircularBuffer{Array{Float64,3}}(data["CircularBufferSize"])
    push!(ADMM["Imbalances"]["EOM"], zeros(data["nTimesteps"], data["nReprDays"], data["nYears"]))
    ADMM["Imbalances"]["H2"] = CircularBuffer{Array{Float64,3}}(data["CircularBufferSize"])
    push!(ADMM["Imbalances"]["H2"], zeros(data["nTimesteps"], data["nReprDays"], data["nYears"]))
    ADMM["Imbalances"]["CM"] = CircularBuffer{Float64}(data["CircularBufferSize"])
    push!(ADMM["Imbalances"]["CM"], 0.0)
    ADMM["Imbalances"]["HCM"] = CircularBuffer{Float64}(data["CircularBufferSize"])
    push!(ADMM["Imbalances"]["HCM"], 0.0)

    ADMM["Residuals"] = Dict()
    ADMM["Residuals"]["Primal"] = Dict()
    ADMM["Residuals"]["Primal"]["EOM"] = CircularBuffer{Float64}(data["CircularBufferSize"])
    push!(ADMM["Residuals"]["Primal"]["EOM"], 0)
    ADMM["Residuals"]["Primal"]["H2"] = CircularBuffer{Float64}(data["CircularBufferSize"])
    push!(ADMM["Residuals"]["Primal"]["H2"], 0)
    ADMM["Residuals"]["Primal"]["CM"] = CircularBuffer{Float64}(data["CircularBufferSize"])
    push!(ADMM["Residuals"]["Primal"]["CM"], 0)
    ADMM["Residuals"]["Primal"]["HCM"] = CircularBuffer{Float64}(data["CircularBufferSize"])
    push!(ADMM["Residuals"]["Primal"]["HCM"], 0)

    ADMM["Residuals"]["Dual"] = Dict()
    ADMM["Residuals"]["Dual"]["EOM"] = CircularBuffer{Float64}(data["CircularBufferSize"])
    push!(ADMM["Residuals"]["Dual"]["EOM"], 0)
    ADMM["Residuals"]["Dual"]["H2"] = CircularBuffer{Float64}(data["CircularBufferSize"])
    push!(ADMM["Residuals"]["Dual"]["H2"], 0)
    ADMM["Residuals"]["Dual"]["CM"] = CircularBuffer{Float64}(data["CircularBufferSize"])
    push!(ADMM["Residuals"]["Dual"]["CM"], 0)
    ADMM["Residuals"]["Dual"]["HCM"] = CircularBuffer{Float64}(data["CircularBufferSize"])
    push!(ADMM["Residuals"]["Dual"]["HCM"], 0)

    ADMM["Tolerance"] = Dict()
    ADMM["Tolerance"]["EOM"] = data["epsilon"] / 1000 * maximum(EOM["D"]) * sqrt(data["nTimesteps"] * data["nReprDays"] * data["nYears"])
    ADMM["Tolerance"]["H2"] = data["epsilon"] / 1000 * maximum(H2["D"]) * sqrt(data["nTimesteps"] * data["nReprDays"] * data["nYears"])
    ADMM["Tolerance"]["CM"] = data["epsilon"] / 1000 * maximum(CM["D"]) * sqrt(1)
    ADMM["Tolerance"]["HCM"] = data["epsilon"] / 1000 * maximum(HCM["D"]) * sqrt(1)

    ADMM["ρ"] = Dict()
    ADMM["ρ"]["EOM"] = CircularBuffer{Float64}(data["CircularBufferSize"])
    push!(ADMM["ρ"]["EOM"], data["rho_EOM"])
    ADMM["ρ"]["H2"] = CircularBuffer{Float64}(data["CircularBufferSize"])
    push!(ADMM["ρ"]["H2"], data["rho_H2"])
    ADMM["ρ"]["CM"] = CircularBuffer{Float64}(data["CircularBufferSize"])
    push!(ADMM["ρ"]["CM"], data["rho_CM"])
    ADMM["ρ"]["HCM"] = CircularBuffer{Float64}(data["CircularBufferSize"])
    push!(ADMM["ρ"]["HCM"], data["rho_HCM"])


    ADMM["n_iter"] = 1
    ADMM["walltime"] = 0

    return results, ADMM
end