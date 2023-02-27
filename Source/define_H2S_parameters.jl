function define_H2S_parameters!(m::String,mod::Model, data::Dict, repr_days::DataFrame, order_matrix::DataFrame)
    # Parameters
    if m == "H2demand"
        mod.ext[:parameters][:WTP] = data["WTP"]
    elseif m == "electrolysis"
        mod.ext[:parameters][:η_E_H2] = data["efficiency_E_H2"] # - 
        mod.ext[:parameters][:IC] = data["OC"]/data["Lifetime"] # MEUR/GW 
    elseif m == "H2storage"
        mod.ext[:sets][:JD_A] = 1:data["nDays"]  # define all periods (days) used for storage
        mod.ext[:sets][:index_repr]= [repr_days[!,:periods][jd] for jd in mod.ext[:sets][:JD]]    # save the day of the year in a compact vector
        mod.ext[:sets][:order_matrix] = Matrix(order_matrix)
        mod.ext[:parameters][:VC] = data["VC"]
        mod.ext[:parameters][:η_ch] = data["efficiency_ch"]
        mod.ext[:parameters][:η_dh] = data["efficiency_dh"]
        mod.ext[:parameters][:IC_cap] = data["OC_cap"]/data["Lifetime"] # MEUR/GW - costs for installed capacity
        mod.ext[:parameters][:IC_vol] = data["OC_vol"]/data["Lifetime"] # MEUR/TWh - costs for installed volume
    end
end