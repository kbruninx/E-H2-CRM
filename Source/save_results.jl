# Save results
function save_results(mdict::Dict,EOM::Dict,H2::Dict,ADMM::Dict,results::Dict,data::Dict,agents::Dict) 
   
    # Power sector
    
    CSV.write(joinpath(home_dir,string("Results_",data["nReprDays"],"_repr_days"),string("demand_ps.csv")), DataFrame(EOM["D"],:auto), delim=";");
    
    cap = zeros(length(agents[:ps]))
    mm = 1
    for m in agents[:ps]
        if m == "Edemand"
            cap[mm] = 0
        else
            cap[mm] = value.(mdict[m].ext[:variables][:cap])
        end
        CSV.write(joinpath(home_dir,string("Results_",data["nReprDays"],"_repr_days"),string("Generation",agents[:ps][mm],".csv")), DataFrame(results["g"][agents[:ps][mm]][end],:auto), delim=";");
        mm = mm + 1
    end
    CSV.write(joinpath(home_dir,string("Results_",data["nReprDays"],"_repr_days"),string("capacity_ps.csv")), DataFrame(transpose(cap),:auto), delim=";",header=string.("CAP_",agents[:ps]));

    ENS = EOM["D"] + results["g"]["Edemand"][end]
    CSV.write(joinpath(home_dir,string("Results_",data["nReprDays"],"_repr_days"),string("ENS_el.csv")), DataFrame(ENS,:auto), delim=";");

    # Electricity prices
    CSV.write(joinpath(home_dir,string("Results_",data["nReprDays"],"_repr_days"),string("electricity_price.csv")), DataFrame(results["λ"]["EOM"][end],:auto), delim=";");

    # Hydrogen prices
    CSV.write(joinpath(home_dir,string("Results_",data["nReprDays"],"_repr_days"),string("hydrogen_price.csv")), DataFrame(results["λ"]["H2"][end],:auto), delim=";");

    # Hydrogen sector
    CSV.write(joinpath(home_dir,string("Results_",data["nReprDays"],"_repr_days"),string("SOC",agents[:h2s][3],".csv")), DataFrame(results["SOC"][end],:auto), delim=";");
    CSV.write(joinpath(home_dir,string("Results_",data["nReprDays"],"_repr_days"),string("SOC_AD_0",agents[:h2s][3],".csv")), DataFrame(results["SOC_AD_0"][end]',:auto), delim=";");
    CSV.write(joinpath(home_dir,string("Results_",data["nReprDays"],"_repr_days"),string("Charge",agents[:h2s][3],".csv")), DataFrame(results["chH"][end],:auto), delim=";");
    CSV.write(joinpath(home_dir,string("Results_",data["nReprDays"],"_repr_days"),string("Discharge",agents[:h2s][3],".csv")), DataFrame(results["dhH"][end],:auto), delim=";");
    CSV.write(joinpath(home_dir,string("Results_",data["nReprDays"],"_repr_days"),string("demand_H2s.csv")), DataFrame(H2["D"],:auto), delim=";");

    mm = 1
    capH = zeros(length(agents[:h2s]))
    volH = zeros(1) # set at one because there is only 1 storage asset, can be incremented
    for m in agents[:h2s]
        if m == "H2demand"
            capH[mm] = 0
        else
            capH[mm] = value.(mdict[m].ext[:variables][:capH])
            if m == "H2storage"
                volH[1] = value.(mdict[m].ext[:variables][:volH])
            end
        end
        CSV.write(joinpath(home_dir,string("Results_",data["nReprDays"],"_repr_days"),string("H2Generation",agents[:h2s][mm],".csv")), DataFrame(results["h2"][agents[:h2s][mm]][end],:auto), delim=";");
        mm = mm + 1  
    end
    CSV.write(joinpath(home_dir,string("Results_",data["nReprDays"],"_repr_days"),string("capacity_h2s.csv")), DataFrame(transpose(capH),:auto), delim=";",header=string.("CAP_",agents[:h2s]));
    CSV.write(joinpath(home_dir,string("Results_",data["nReprDays"],"_repr_days"),string("storage_volume.csv")), DataFrame(transpose(volH),:auto), delim=";");

end
