# Save results
function save_results(mdict::Dict,EOM::Dict,H2::Dict,ADMM::Dict,results::Dict,data::Dict,agents::Dict) 
    # Power sector
    cap = zeros(length(agents[:ps]))
    mm = 1
    for m in agents[:ps]
        cap[mm] = value.(mdict[m].ext[:variables][:cap])
        mm = mm+1
    end
    CSV.write(joinpath(home_dir,string("Results_",data["nReprDays"],"_repr_days"),string("capacity_ps.csv")), DataFrame(transpose(cap),:auto), delim=";",header=string.("CAP_",agents[:ps]));

    # Electricity prices
    CSV.write(joinpath(home_dir,string("Results_",data["nReprDays"],"_repr_days"),string("electricity_price.csv")), DataFrame(results["λ"]["EOM"][end],:auto), delim=";");

    # Hydrogen sector
    CSV.write(joinpath(home_dir,string("Results_",data["nReprDays"],"_repr_days"),string("hydrogen_production_",agents[:h2s][1],".csv")), DataFrame(results["h2"][agents[:h2s][1]][end],:auto), delim=";");
    CSV.write(joinpath(home_dir,string("Results_",data["nReprDays"],"_repr_days"),string("hydrogen_exchange_",agents[:h2s][2],".csv")), DataFrame(results["h2"][agents[:h2s][2]][end],:auto), delim=";");
    
    mm = 1
    capH = zeros(length(agents[:h2s]))
    for m in agents[:h2s]
        capH[mm] = value.(mdict[m].ext[:variables][:capH])
        mm = mm + 1
    end
    CSV.write(joinpath(home_dir,string("Results_",data["nReprDays"],"_repr_days"),string("capacity_h2s.csv")), DataFrame(transpose(capH),:auto), delim=";",header=string.("CAP_",agents[:h2s]));

end
