function ADMM_subroutine!(m::String,data::Dict,results::Dict,ADMM::Dict,EOM::Dict,H2::Dict,mod::Model,agents::Dict,TO::TimerOutput)
TO_local = TimerOutput()
# Calculate penalty terms ADMM and update price to most recent value 
@timeit TO_local "Compute ADMM penalty terms" begin
    if mod.ext[:parameters][:EOM] == 1
        mod.ext[:parameters][:g_bar] = results["g"][m][end] - 1/(EOM["nAgents"]+1)*ADMM["Imbalances"]["EOM"][end]
        mod.ext[:parameters][:λ_EOM] = results["λ"]["EOM"][end] 
        mod.ext[:parameters][:ρ_EOM] = ADMM["ρ"]["EOM"][end]
    end
    if mod.ext[:parameters][:H2] == 1
        mod.ext[:parameters][:gH_bar] = results["h2"][m][end] - 1/(H2["nAgents"]+1)*ADMM["Imbalances"]["H2"][end]
        mod.ext[:parameters][:λ_H2] = results["λ"]["H2"][end] 
        mod.ext[:parameters][:ρ_H2] = ADMM["ρ"]["H2"][end]
    end
end

# Solve agents decision problems:
if m in agents[:ps]
    @timeit TO_local "Solve power sector" begin
        solve_ps_agent!(m,mod)  
    end
elseif m in agents[:h2s]
    @timeit TO_local "Solve hydrogen sector" begin
        solve_h2s_agent!(m,mod)  
    end
end

# Query results
@timeit TO_local "Query results" begin
    if mod.ext[:parameters][:EOM] == 1
        push!(results["g"][m], collect(value.(mod.ext[:variables][:g])))
    end
    if mod.ext[:parameters][:H2] == 1
        push!(results["h2"][m], collect(value.(mod.ext[:expressions][:gH])))
    end                     
end

# Merge local TO with TO:
merge!(TO,TO_local)
end