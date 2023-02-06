function define_results!(data::Dict, results::Dict, ADMM::Dict, agents::Dict, EOM::Dict, H2::Dict)
    results["g"] = Dict()
    results["h2"] = Dict()
    results["e_demand"] = Dict()
    results["h2_demand"] = Dict()
    results["SOC"] = Dict()

    for m in agents[:eom]
        results["g"][m] = CircularBuffer{Array{Float64,2}}(data["CircularBufferSize"])
        push!(results["g"][m], zeros(data["nTimesteps"], data["nReprDays"]))
    end
    for m in agents[:h2]
        results["h2"][m] = CircularBuffer{Array{Float64,2}}(data["CircularBufferSize"])
        push!(results["h2"][m], zeros(data["nTimesteps"], data["nReprDays"]))
        if m == "H2storage"
            results["SOC"] = CircularBuffer{Array{Float64,2}}(data["CircularBufferSize"])
            push!(results["SOC"], zeros(data["nTimesteps"], data["nReprDays"]))
            results["dhH"] = CircularBuffer{Array{Float64,2}}(data["CircularBufferSize"])
            push!(results["dhH"], zeros(data["nTimesteps"], data["nReprDays"]))
            results["chH"] = CircularBuffer{Array{Float64,2}}(data["CircularBufferSize"])
            push!(results["chH"], zeros(data["nTimesteps"], data["nReprDays"]))
        end
    end

    results["λ"] = Dict()
    results["λ"]["EOM"] = CircularBuffer{Array{Float64,2}}(data["CircularBufferSize"])
    push!(results["λ"]["EOM"], zeros(data["nTimesteps"], data["nReprDays"]))
    results["λ"]["H2"] = CircularBuffer{Array{Float64,2}}(data["CircularBufferSize"])
    push!(results["λ"]["H2"], zeros(data["nTimesteps"], data["nReprDays"]))

    ADMM["Imbalances"] = Dict()
    ADMM["Imbalances"]["EOM"] = CircularBuffer{Array{Float64,2}}(data["CircularBufferSize"])
    push!(ADMM["Imbalances"]["EOM"], zeros(data["nTimesteps"], data["nReprDays"]))
    ADMM["Imbalances"]["H2"] = CircularBuffer{Array{Float64,2}}(data["CircularBufferSize"])
    push!(ADMM["Imbalances"]["H2"], zeros(data["nTimesteps"], data["nReprDays"]))

    ADMM["Residuals"] = Dict()
    ADMM["Residuals"]["Primal"] = Dict()
    ADMM["Residuals"]["Primal"]["EOM"] = CircularBuffer{Float64}(data["CircularBufferSize"])
    push!(ADMM["Residuals"]["Primal"]["EOM"], 0)
    ADMM["Residuals"]["Primal"]["H2"] = CircularBuffer{Float64}(data["CircularBufferSize"])
    push!(ADMM["Residuals"]["Primal"]["H2"], 0)

    ADMM["Residuals"]["Dual"] = Dict()
    ADMM["Residuals"]["Dual"]["EOM"] = CircularBuffer{Float64}(data["CircularBufferSize"])
    push!(ADMM["Residuals"]["Dual"]["EOM"], 0)
    ADMM["Residuals"]["Dual"]["H2"] = CircularBuffer{Float64}(data["CircularBufferSize"])
    push!(ADMM["Residuals"]["Dual"]["H2"], 0)

    ADMM["Tolerance"] = Dict()
    ADMM["Tolerance"]["EOM"] = data["epsilon"] / 100 * maximum(EOM["D"]) * sqrt(data["nTimesteps"] * data["nReprDays"])
    ADMM["Tolerance"]["H2"] = data["epsilon"] / 100 * maximum(H2["D"]) * sqrt(data["nTimesteps"] * data["nReprDays"])

    ADMM["ρ"] = Dict()
    ADMM["ρ"]["EOM"] = CircularBuffer{Float64}(data["CircularBufferSize"])
    push!(ADMM["ρ"]["EOM"], data["rho_EOM"])
    ADMM["ρ"]["H2"] = CircularBuffer{Float64}(data["CircularBufferSize"])
    push!(ADMM["ρ"]["H2"], data["rho_H2"])

    ADMM["n_iter"] = 1
    ADMM["walltime"] = 0

    return results, ADMM
end