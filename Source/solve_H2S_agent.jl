function solve_h2s_agent!(m::String, mod::Model, H2::Dict)
   # Extract sets
   JH = mod.ext[:sets][:JH]
   JD = mod.ext[:sets][:JD]
   JY = mod.ext[:sets][:JY]

   # Extract common parameters
   W = mod.ext[:parameters][:W] # weight of the representative days
   P = mod.ext[:parameters][:P] # probability of each scenario
   γ = mod.ext[:parameters][:γ] # weight of expected revenues and CVAR
   β = mod.ext[:parameters][:β] # risk aversion parametrization

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
      α = mod.ext[:variables][:α]
      u = mod.ext[:variables][:u]

      # Expressions

      gH_positive = mod.ext[:expressions][:gH_positive] = @expression(mod, gH_VOLL + gH_ela)
      gH = mod.ext[:expressions][:gH] = @expression(mod, - gH_positive)

      profit = mod.ext[:expressions][:profit] = @expression(mod, [jy = JY],      # the welfare is formulated negative to make the problem convex
      sum(W[jd, jy] * (WTP * gH_positive[jh, jd, jy] - (gH_ela[jh, jd, jy])^2 * WTP / (2 * ela_H2[jy])
      - λ_H2[jh, jd, jy] * gH_positive[jh, jd, jy]) for jh in JH, jd in JD) )

      CVAR = mod.ext[:expressions][:CVAR] = @expression(mod,
      α - ((1/β) * sum(P[jy] * u[jy] for jy in JY)) )

      # Update objective function 
      mod.ext[:objective] = @objective(mod, Min,  
      - γ * sum(P[jy] * profit[jy] for jy in JY) - (1 - γ) * CVAR
      + sum(ρ_H2 / 2 * W[jd, jy] * (gH[jh, jd, jy] - gH_bar[jh, jd, jy])^2 for jh in JH, jd in JD, jy in JY)
      ) 

      if γ < 1
         for jy in JY
            delete(mod, mod.ext[:constraints][:VAR_threshold][jy])
         end
         
         # CVAR constraint
         mod.ext[:constraints][:VAR_threshold] = @constraint(mod, [jy = JY],
         α - profit[jy] <= u[jy] )
      end

   elseif m == "electrolysis"

      # Extract parameters
      IC = mod.ext[:parameters][:IC] # annuity investment costs
      η_E_H2 = mod.ext[:parameters][:η_E_H2] # efficiency electrolysis

      # Decision variables
      capH = mod.ext[:variables][:capH]
      g = mod.ext[:variables][:g]  # note this is defined as a negative number, consumption
      α = mod.ext[:variables][:α]
      u = mod.ext[:variables][:u]

      # Create affine expressions  
      #inv_cost = mod.ext[:expressions][:inv_cost]
      
      gH = mod.ext[:expressions][:gH] = @expression(mod, -η_E_H2 * g) # [TWh]
      
      cost = mod.ext[:expressions][:cost] = @expression(mod, [jy = JY],
      IC * capH
      - sum(W[jd, jy] * λ_EOM[jh, jd, jy] * g[jh, jd, jy] for jh in JH, jd in JD))
      
      revenue = mod.ext[:expressions][:revenue] = @expression(mod, [jy = JY],
      sum(W[jd, jy] * λ_H2[jh, jd, jy] * gH[jh, jd, jy] for jh in JH, jd in JD))

      profit = mod.ext[:expressions][:profit] = @expression(mod, [jy = JY],
      revenue[jy] - cost[jy] )

      CVAR = mod.ext[:expressions][:CVAR] = @expression(mod,
      α - ((1/β) * sum(P[jy] * u[jy] for jy in JY)) )

      #tot_revenue = mod.ext[:expressions][:tot_revenue] = @expression(mod,
      #sum(P[jy] * sum(W[jd, jy] * λ_H2[jh, jd, jy] * gH[jh, jd, jy] for jh in JH, jd in JD) for jy in JY))

      # Update objective function
      mod.ext[:objective] = @objective(mod, Min,
      - γ * sum(P[jy] * (revenue[jy] - cost[jy]) for jy in JY) - (1 - γ) * CVAR
      + sum(ρ_EOM / 2 * W[jd, jy] * (g[jh, jd, jy] - g_bar[jh, jd, jy])^2 for jh in JH, jd in JD, jy in JY)
      + sum(ρ_H2 / 2 * W[jd, jy] * (gH[jh, jd, jy] - gH_bar[jh, jd, jy])^2 for jh in JH, jd in JD, jy in JY)
      )

      if γ < 1
         for jy in JY
            delete(mod, mod.ext[:constraints][:VAR_threshold][jy])
         end

         # CVAR constraint
         mod.ext[:constraints][:VAR_threshold] = @constraint(mod, [jy = JY],
         α - (revenue[jy] - cost[jy]) <= u[jy] )
      end      

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
      α = mod.ext[:variables][:α]
      u = mod.ext[:variables][:u]

      # Create affine expressions  
      #inv_cost = mod.ext[:expressions][:inv_cost]
      
      gH = mod.ext[:expressions][:gH] = @expression(mod, dhH - chH)
      
      cost = mod.ext[:expressions][:cost] = @expression(mod, [jy = JY],
      IC_cap * capH + IC_vol * volH
      + sum(W[jd, jy] * λ_H2[jh, jd, jy] * chH[jh, jd, jy] for jh in JH, jd in JD))
      
      revenue = mod.ext[:expressions][:revenue] = @expression(mod, [jy = JY],
      sum(W[jd, jy] * λ_H2[jh, jd, jy] * dhH[jh, jd, jy] for jh in JH, jd in JD))

      profit = mod.ext[:expressions][:profit] = @expression(mod, [jy = JY],
      revenue[jy] - cost[jy] )

      CVAR = mod.ext[:expressions][:CVAR] = @expression(mod,
      α - ((1/β) * sum(P[jy] * u[jy] for jy in JY)) )

      #tot_revenue = mod.ext[:expressions][:tot_revenue] = @expression(mod,
      #sum(P[jy] * sum(W[jd, jy] * λ_H2[jh, jd, jy] * dhH[jh, jd, jy] for jh in JH, jd in JD) for jy in JY))


      # Update objective function
      mod.ext[:objective] = @objective(mod, Min,
      - γ * sum(P[jy] * (revenue[jy] - cost[jy]) for jy in JY) - (1 - γ) * CVAR
      + sum(ρ_H2 / 2 * W[jd, jy] * (gH[jh, jd, jy] - gH_bar[jh, jd, jy])^2 for jh in JH, jd in JD, jy in JY)
      )

      if γ < 1
         for jy in JY
            delete(mod, mod.ext[:constraints][:VAR_threshold][jy])
         end

         # CVAR constraint
         mod.ext[:constraints][:VAR_threshold] = @constraint(mod, [jy = JY],
         α - (revenue[jy] - cost[jy]) <= u[jy] )
      end
   end

   # solve problem
   optimize!(mod)

   return mod

end