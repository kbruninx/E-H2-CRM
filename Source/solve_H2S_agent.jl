function solve_h2s_agent!(m::String,mod::Model,H2::Dict)
   # Extract sets
   JH = mod.ext[:sets][:JH]
   JD = mod.ext[:sets][:JD]

   # Extract common parameters
   W = mod.ext[:parameters][:W] # weight of the representative days

   # ADMM algorithm parameters
   λ_EOM = mod.ext[:parameters][:λ_EOM] # EOM prices
   g_bar = mod.ext[:parameters][:g_bar] # element in ADMM penalty term related to EOM
   ρ_EOM = mod.ext[:parameters][:ρ_EOM] # rho-value in ADMM related to EUA auctions
   λ_H2 = mod.ext[:parameters][:λ_H2] # H2 prices
   gH_bar = mod.ext[:parameters][:gH_bar] # element in ADMM penalty term related to hydrogen market
   ρ_H2 = mod.ext[:parameters][:ρ_H2] # rho-value in ADMM related to H2 market

   if m == "H2demand"
      # Extract parameters
      WTP = mod.ext[:parameters][:WTP]  # willingness to pay ("price cap") of consumers
      ela_H2 = H2["elasticity"] # section of price-elasticity [MWh]

      # Decision variables/expressions - all negative here
      gH_VOLL = mod.ext[:variables][:gH_VOLL]
      gH_ela = mod.ext[:variables][:gH_ela]
      gH = mod.ext[:expressions][:gH]

      # Update objective function 
      mod.ext[:objective] = @objective(mod, Min,   # minimize a negative quantity (maximize its absolute value)
      sum(W[jd] * ((WTP * gH[jh,jd] + (gH_ela[jh,jd])^2 * WTP / (2*ela_H2)) - λ_H2[jh, jd] * gH[jh,jd]) for jh in JH, jd in JD)
      + sum(ρ_H2 / 2 * W[jd] * (gH[jh, jd] - gH_bar[jh, jd])^2 for jh in JH, jd in JD)
      )
      # gH_ela is defined negative, therefore the function here must be adapted since there is gH_ela^2 : put + instead of - in front of gH_ela^2

   elseif m == "electrolysis"  # NOT SURE IT WORKS

      # Extract parameters
      IC = mod.ext[:parameters][:IC] # annuity investment costs
      η_E_H2 = mod.ext[:parameters][:η_E_H2] # efficiency electrolysis

      # Decision variables
      capH = mod.ext[:variables][:capH] 
      g = mod.ext[:variables][:g]  # note this is defined as a negative number, consumption

      # Create affine expressions  
      inv_cost = mod.ext[:expressions][:inv_cost]
      gH = mod.ext[:expressions][:gH] # [TWh]

      # Update objective function
      mod.ext[:objective] = @objective(mod, Min,
         + inv_cost # [MEUR]
         - sum(W[jd] * (λ_EOM[jh, jd]) * g[jh, jd] for jh in JH, jd in JD) # [MEUR]
         - sum(W[jd] * λ_H2[jh, jd] * gH[jh, jd] for jh in JH, jd in JD)
         + sum(ρ_EOM / 2 * W[jd] * (g[jh, jd] - g_bar[jh, jd])^2 for jh in JH, jd in JD)
         + sum(ρ_H2 / 2 * W[jd] * (gH[jh, jd] - gH_bar[jh, jd])^2 for jh in JH, jd in JD)
      )

   elseif m == "H2storage"

      # Extract parameters
      VC = mod.ext[:parameters][:VC] # variable cost storage
      η_ch = mod.ext[:parameters][:η_ch] # charging efficiency 
      η_dh = mod.ext[:parameters][:η_dh] # discharging efficiency
      IC_cap = mod.ext[:parameters][:IC_cap] # annuity investment cost for capacity
      IC_vol = mod.ext[:parameters][:IC_vol] # annuity investment cost for volume

      # Decision variables
      capH = mod.ext[:variables][:capH] 
      volH = mod.ext[:variables][:volH] 
      chH = mod.ext[:variables][:chH] 
      dhH = mod.ext[:variables][:dhH] 
      SOC = mod.ext[:variables][:SOC]

      # Create affine expressions  
      inv_cost = mod.ext[:expressions][:inv_cost]
      # SOC = mod.ext[:expressions][:SOC] = @expression(mod, SOC + η_ch * chH - dhH / η_dh) # state of charge update
      gH = mod.ext[:expressions][:gH] 

      # Update objective function
      mod.ext[:objective] = @objective(mod, Min,
         + inv_cost # [MEUR]
         - sum(W[jd] * (λ_H2[jh, jd]) * dhH[jh, jd] for jh in JH, jd in JD) # [MEUR]
         + sum(W[jd] * λ_H2[jh, jd] * chH[jh, jd] for jh in JH, jd in JD)
         + sum(ρ_H2 / 2 * W[jd] * (gH[jh, jd] - gH_bar[jh, jd])^2 for jh in JH, jd in JD)
      )
   end
    
   # solve problem
   optimize!(mod);

   return mod

end