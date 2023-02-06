function define_H2S_parameters!(m::String,mod::Model, data::Dict)
    # Parameters 
    mod.ext[:parameters][:IC] = data["OC"]/data["Lifetime"] # EUR/MW or MEUR/TW
    
    if m == "electrolysis"
        mod.ext[:parameters][:η_E_H2] = data["efficiency_E_H2"] # - 
    elseif m == "H2storage"
        mod.ext[:parameters][:VC] = data["VC"]
        mod.ext[:parameters][:η_ch] = data["efficiency_ch"]
        mod.ext[:parameters][:η_dh] = data["efficiency_dh"]
    #elseif m == "Edemand"
     #   mod.ext[:parameters][:WTP] = data["WTP"]
    end
end