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

    # Hydrogen production
    CSV.write(joinpath(home_dir,string("Results_",data["nReprDays"],"_repr_days"),string("hydrogen_production_",agents[:h2s][1],".csv")), DataFrame(results["h2"][agents[:h2s][1]][end],:auto), delim=";");
end
