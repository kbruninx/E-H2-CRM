function update_rho!(ADMM::Dict, iter::Int64)
    if mod(iter,1) == 0 # can be used to do update ever i-th iteration
        # ρ-updates following Boyd et al.  
        if ADMM["Residuals"]["Primal"]["EOM"][end] > 2*ADMM["Residuals"]["Dual"]["EOM"][end]
            push!(ADMM["ρ"]["EOM"], minimum([1000,1.1*ADMM["ρ"]["EOM"][end]]))
        elseif ADMM["Residuals"]["Dual"]["EOM"][end] > 2*ADMM["Residuals"]["Primal"]["EOM"][end]
            push!(ADMM["ρ"]["EOM"], 1/1.1*ADMM["ρ"]["EOM"][end])
        end

        if ADMM["Residuals"]["Primal"]["H2"][end] > 2*ADMM["Residuals"]["Dual"]["H2"][end]
            push!(ADMM["ρ"]["H2"], minimum([1000,1.1*ADMM["ρ"]["H2"][end]]))
        elseif ADMM["Residuals"]["Dual"]["H2"][end] > 2*ADMM["Residuals"]["Primal"]["H2"][end]
            push!(ADMM["ρ"]["H2"], 1/1.1*ADMM["ρ"]["H2"][end])
        end
    end
end