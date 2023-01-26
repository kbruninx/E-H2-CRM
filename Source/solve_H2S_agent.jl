function solve_h2s_agent!(mod::Model)
   # Extract sets
   JH = mod.ext[:sets][:JH]
   JD = mod.ext[:sets][:JD]
      
   # Extract parameters 
   W = mod.ext[:parameters][:W] # weight of the representative days
   IC = mod.ext[:parameters][:IC] # overnight investment costs
   λ_EOM = mod.ext[:parameters][:λ_EOM] # EOM prices
   ρ_EOM = mod.ext[:parameters][:ρ_EOM] # rho-value in ADMM related to EUA auctions
   g_bar = mod.ext[:parameters][:g_bar] # element in ADMM penalty term related to EOM
   λ_H2 = mod.ext[:parameters][:λ_H2] # H2 prices
   ρ_H2 = mod.ext[:parameters][:ρ_H2] # rho-value in ADMM related to H2 market
   gH_bar = mod.ext[:parameters][:gH_bar] # element in ADMM penalty term related to hydrogen market

   # Extract variables and expressions
   capH = mod.ext[:variables][:capH] 
   g = mod.ext[:variables][:g] 
   gH = mod.ext[:variables][:gH]

   # Update objective
   mod.ext[:objective] = @objective(mod, Min,
        + IC*capH # [MEUR]
        - sum(W[jd]*(λ_EOM[jh,jd])*g[jh,jd] for jh in JH, jd in JD) # [MEUR]
        - sum(W[jd]*λ_H2[jh,jd]*gH[jh,jd] for jh in JH, jd in JD)
        + sum(ρ_EOM/2*W[jd]*(g[jh,jd] - g_bar[jh,jd])^2 for jh in JH, jd in JD)  
        + sum(ρ_H2/2*W[jd]*(gH[jh,jd] - gH_bar[jh,jd])^2 for jh in JH, jd in JD) 
   )
    
   # solve problem
   optimize!(mod);

   return mod

end