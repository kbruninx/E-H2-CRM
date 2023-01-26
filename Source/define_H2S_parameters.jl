function define_H2S_parameters!(mod::Model, data::Dict)
    # Parameters 
    mod.ext[:parameters][:η_E_H2] = data["efficiency_E_H2"] # - 
    mod.ext[:parameters][:IC] = data["OC"]/data["Lifetime"] # EUR/MW or MEUR/TW
    
    return mod
end