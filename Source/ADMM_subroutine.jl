function ADMM_subroutine!(m::String,data::Dict,results::Dict,ADMM::Dict,EOM::Dict,H2::Dict,mod::Model,agents::Dict,TO::TimerOutput)
TO_local = TimerOutput()
# Calculate penalty terms ADMM and update price to most recent value 
@timeit TO_local "Compute ADMM penalty terms" begin
    if mod.ext[:parameters][:EOM] == 1
        mod.ext[:parameters][:g_bar] = results["g"][m][end] - 1 / (EOM["nAgents"] + 1) * ADMM["Imbalances"]["EOM"][end]
        mod.ext[:parameters][:λ_EOM] = results["λ"]["EOM"][end] 
        mod.ext[:parameters][:ρ_EOM] = ADMM["ρ"]["EOM"][end]
    end
    if mod.ext[:parameters][:H2] == 1
        mod.ext[:parameters][:gH_bar] = results["h2"][m][end] - 1/(H2["nAgents"]+1)*ADMM["Imbalances"]["H2"][end]
        mod.ext[:parameters][:λ_H2] = results["λ"]["H2"][end] 
        mod.ext[:parameters][:ρ_H2] = ADMM["ρ"]["H2"][end]
    end
    if mod.ext[:parameters][:CM] == 1
        mod.ext[:parameters][:cap_bar] = results["cap_cm"][m][end] - 1/(CM["nAgents"]+1)*ADMM["Imbalances"]["CM"][end]
        mod.ext[:parameters][:λ_CM] = results["λ"]["CM"][end] 
        mod.ext[:parameters][:ρ_CM] = ADMM["ρ"]["CM"][end]
    end
    if mod.ext[:parameters][:HCM] == 1
        mod.ext[:parameters][:capH_bar] = results["capH_cm"][m][end] - 1 / (HCM["nAgents"] + 1) * ADMM["Imbalances"]["HCM"][end]
        mod.ext[:parameters][:λ_HCM] = results["λ"]["HCM"][end]
        mod.ext[:parameters][:ρ_HCM] = ADMM["ρ"]["HCM"][end]
    end
end

# Solve agents decision problems:
if m in agents[:ps]
    @timeit TO_local "Solve power sector" begin
        solve_ps_agent!(m,mod,EOM)  
    end
elseif m in agents[:h2s]
    @timeit TO_local "Solve hydrogen sector" begin
        solve_h2s_agent!(m,mod,H2)  
    end
end

# Query results
@timeit TO_local "Query results" begin
    if mod.ext[:parameters][:EOM] == 1
        if m == "Edemand"
            push!(results["g"][m], collect(value.(mod.ext[:expressions][:g])))
            push!(results["g_ela"], collect(value.(mod.ext[:variables][:g_ela])))
        else
            push!(results["g"][m], collect(value.(mod.ext[:variables][:g])))
        end
    end
    if mod.ext[:parameters][:H2] == 1
        push!(results["h2"][m], collect(value.(mod.ext[:expressions][:gH])))
        if m == "H2storage"
            push!(results["SOC"], collect(value.(mod.ext[:variables][:SOC])))
            push!(results["SOC_AD_0"], collect(value.(mod.ext[:variables][:SOC_AD_0])))
            push!(results["dhH"], collect(value.(mod.ext[:variables][:dhH])))
            push!(results["chH"], collect(value.(mod.ext[:variables][:chH])))
        elseif m == "H2demand"
            push!(results["gH2_ela"], collect(value.(mod.ext[:variables][:gH_ela])))
        end
    end     
    if mod.ext[:parameters][:CM] == 1
        push!(results["cap_cm"][m], first(collect(value.(mod.ext[:variables][:cap_cm]))))
    end
    if mod.ext[:parameters][:HCM] == 1
        push!(results["capH_cm"][m], first(collect(value.(mod.ext[:variables][:capH_cm]))))
    end
end

# Merge local TO with TO:
merge!(TO,TO_local)
end