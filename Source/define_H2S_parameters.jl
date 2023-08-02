function define_H2S_parameters!(m::String, mod::Model, data::Dict, repr_days::Dict, order_matrix::Dict)
    # Parameters
    if m == "H2demand"
        mod.ext[:parameters][:WTP] = data["WTP"]
    elseif m == "electrolysis"
        mod.ext[:parameters][:η_E_H2] = data["efficiency_E_H2"] # - 
        mod.ext[:parameters][:IC] = (data["discount_rate"] * data["OC"]) / (1 - (1 / ((1 + data["discount_rate"])^data["Lifetime"]))) # MEUR/GW 
    elseif m == "H2storage"
        mod.ext[:sets][:JD_A] = 1:data["nDays"]  # define all periods (days) used for storage
        mod.ext[:parameters][:VC] = data["VC"]
        mod.ext[:parameters][:η_ch] = data["efficiency_ch"]
        mod.ext[:parameters][:η_dh] = data["efficiency_dh"]
        mod.ext[:parameters][:IC_cap] = (data["discount_rate"] * data["OC_cap"]) / (1 - (1 / ((1 + data["discount_rate"])^data["Lifetime"]))) # MEUR/GW - costs for installed capacity
        mod.ext[:parameters][:IC_vol] = (data["discount_rate"] * data["OC_vol"]) / (1 - (1 / ((1 + data["discount_rate"])^data["Lifetime"]))) # MEUR/TWh - costs for installed volume

        mod.ext[:sets][:index_repr]= [repr_days[jy][!,:periods][jd] for jd in mod.ext[:sets][:JD], jy in mod.ext[:sets][:JY]]    # save the day of the year in a compact vector
        
        mod.ext[:sets][:order_matrix] = zeros(size(Matrix(order_matrix[1]),1), size(Matrix(order_matrix[1]),2), data["nYears"])   # initailizing the 3D matrix for days VS repr_days for the different years
        for jy in mod.ext[:sets][:JY]    
            mod.ext[:sets][:order_matrix][:,:,jy] = Matrix(order_matrix[jy])
        end
    end
end