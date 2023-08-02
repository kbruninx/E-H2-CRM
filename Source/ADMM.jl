# ADMM 
function ADMM!(results::Dict,ADMM::Dict,EOM::Dict,H2::Dict,mdict::Dict,agents::Dict,data::Dict,TO::TimerOutput)
    convergence = 0
    σ = mdict["WindOnshore"].ext[:parameters][:σ] # switch capacity market
    σH = mdict["WindOnshore"].ext[:parameters][:σH] # switch capacity market 
    iterations = ProgressBar(1:data["ADMM"]["max_iter"])
    for iter in iterations
        if convergence == 0
            # Multi-threaded version
            @sync for m in agents[:all] 
                # created subroutine to allow multi-treading to solve agents' decision problems
                @spawn ADMM_subroutine!(m,data,results,ADMM,EOM,H2,mdict[m],agents,TO)
            end

            # Single-threaded version
            # for m in agents[:all] 
            #     # created subroutine to allow multi-treading to solve agents' decision problems
            #     ADMM_subroutine!(m,results,ADMM,ETS,EOM,REC,H2,H2CN_prod,H2CN_cap,NG,mdict[m],agents,TO)
            # end

            # Imbalances 
            @timeit TO "Compute imbalances" begin
                push!(ADMM["Imbalances"]["EOM"], sum(results["g"][m][end] for m in agents[:eom]))
                push!(ADMM["Imbalances"]["H2"], sum(results["h2"][m][end] for m in agents[:h2]))
                push!(ADMM["Imbalances"]["CM"], sum(results["cap_cm"][m][end] for m in agents[:cm]) - CM["D"])
                push!(ADMM["Imbalances"]["HCM"], sum(results["capH_cm"][m][end] for m in agents[:hcm]) - HCM["D"])
            end

            # Primal residuals 
            @timeit TO "Compute primal residuals" begin
                push!(ADMM["Residuals"]["Primal"]["EOM"], sqrt(sum(ADMM["Imbalances"]["EOM"][end].^2)))
                push!(ADMM["Residuals"]["Primal"]["H2"], sqrt(sum(ADMM["Imbalances"]["H2"][end].^2)))
                push!(ADMM["Residuals"]["Primal"]["CM"], sqrt(sum(ADMM["Imbalances"]["CM"][end].^2)))
                push!(ADMM["Residuals"]["Primal"]["HCM"], sqrt(sum(ADMM["Imbalances"]["HCM"][end] .^ 2)))
            end

            # Dual residuals
            @timeit TO "Compute dual residuals" begin 
            if iter > 1
                push!(ADMM["Residuals"]["Dual"]["EOM"], sqrt(sum(sum((ADMM["ρ"]["EOM"][end]*((results["g"][m][end]-sum(results["g"][mstar][end] for mstar in agents[:eom])./(EOM["nAgents"]+1)) - (results["g"][m][end-1]-sum(results["g"][mstar][end-1] for mstar in agents[:eom])./(EOM["nAgents"]+1)))).^2 for m in agents[:eom]))))               
                push!(ADMM["Residuals"]["Dual"]["H2"], sqrt(sum(sum((ADMM["ρ"]["H2"][end]*((results["h2"][m][end]-sum(results["h2"][mstar][end] for mstar in agents[:h2])./(H2["nAgents"]+1)) - (results["h2"][m][end-1]-sum(results["h2"][mstar][end-1] for mstar in agents[:h2])./(H2["nAgents"]+1)))).^2 for m in agents[:h2]))))
                push!(ADMM["Residuals"]["Dual"]["CM"], sqrt(sum(sum((ADMM["ρ"]["CM"][end]*((results["cap_cm"][m][end]-sum(results["cap_cm"][mstar][end] for mstar in agents[:cm])./(CM["nAgents"]+1)) - (results["cap_cm"][m][end-1]-sum(results["cap_cm"][mstar][end-1] for mstar in agents[:cm])./(CM["nAgents"]+1)))).^2 for m in agents[:cm]))))
                push!(ADMM["Residuals"]["Dual"]["HCM"], sqrt(sum(sum((ADMM["ρ"]["HCM"][end] * ((results["capH_cm"][m][end] - sum(results["capH_cm"][mstar][end] for mstar in agents[:hcm]) ./ (HCM["nAgents"] + 1)) - (results["capH_cm"][m][end-1] - sum(results["capH_cm"][mstar][end-1] for mstar in agents[:hcm]) ./ (HCM["nAgents"] + 1)))) .^ 2 for m in agents[:hcm]))))
            end
            end

            # Price updates 
            # In general, price caps or floors can be imposed, but may slow down convergence (a negative price gives a stronger incentive than a zero price).
            @timeit TO "Update prices" begin
                push!(results["λ"]["EOM"], results[ "λ"]["EOM"][end] - ADMM["ρ"]["EOM"][end]*ADMM["Imbalances"]["EOM"][end])
                push!(results["λ"]["H2"], results[ "λ"]["H2"][end] - ADMM["ρ"]["H2"][end]*ADMM["Imbalances"]["H2"][end])
                push!(results["λ"]["CM"], results[ "λ"]["CM"][end] - ADMM["ρ"]["CM"][end]*ADMM["Imbalances"]["CM"][end])
                push!(results["λ"]["HCM"], results["λ"]["HCM"][end] - ADMM["ρ"]["HCM"][end] * ADMM["Imbalances"]["HCM"][end])
            end

            # Update ρ-values
            @timeit TO "Update ρ" begin
                update_rho!(ADMM,iter)
            end

            # Progress bar
            @timeit TO "Progress bar" begin
                set_description(iterations, string(@sprintf("ΔEOM %.3f -- ΔH2-h %.3f -- ΔCM %.3f -- ΔHCM %.3f ", ADMM["Residuals"]["Primal"]["EOM"][end], ADMM["Residuals"]["Primal"]["H2"][end], σ * ADMM["Residuals"]["Primal"]["CM"][end], σH * ADMM["Residuals"]["Primal"]["HCM"][end])))
            end

            # Check convergence: primal and dual satisfy tolerance 
            if ADMM["Residuals"]["Primal"]["EOM"][end] <= ADMM["Tolerance"]["EOM"] && ADMM["Residuals"]["Dual"]["EOM"][end] <= ADMM["Tolerance"]["EOM"] && ADMM["Residuals"]["Primal"]["H2"][end] <= ADMM["Tolerance"]["H2"] && ADMM["Residuals"]["Dual"]["H2"][end] <= ADMM["Tolerance"]["H2"] && (σ * ADMM["Residuals"]["Primal"]["CM"][end]) <= ADMM["Tolerance"]["CM"] && (σ * ADMM["Residuals"]["Dual"]["CM"][end]) <= ADMM["Tolerance"]["CM"] && (σH * ADMM["Residuals"]["Primal"]["HCM"][end]) <= ADMM["Tolerance"]["HCM"] && (σH * ADMM["Residuals"]["Dual"]["HCM"][end]) <= ADMM["Tolerance"]["HCM"]
                convergence = 1
            end

            # store number of iterations
            ADMM["n_iter"] = copy(iter)
        end
    end
end
