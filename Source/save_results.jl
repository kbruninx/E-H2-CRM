# Save results
function save_results(mdict::Dict,EOM::Dict,H2::Dict,ADMM::Dict,results::Dict,data::Dict,agents::Dict) 
   
    # Power sector
    for jy in 1:data["nYears"]
        CSV.write(joinpath(home_dir,string("Results_",data["nReprDays"],"_repr_days"),string("demand_ps_",jy,".csv")), DataFrame(EOM["D"][:,:,jy],:auto), delim=";");
    end

    cap = zeros(length(agents[:ps]))
    tot_revenue = zeros(length(agents[:ps]))
    CVAR = zeros(length(agents[:ps]))
    VAR = zeros(length(agents[:ps]))
    mm = 1
    for m in agents[:ps]

        CSV.write(joinpath(home_dir,string("Results_",data["nReprDays"],"_repr_days"),string("Diff_VAR",agents[:ps][mm],".csv")),DataFrame(transpose(value.(mdict[m].ext[:variables][:u]).data),:auto), delim=";");
        CSV.write(joinpath(home_dir,string("Results_",data["nReprDays"],"_repr_days"),string("Profit",agents[:ps][mm],".csv")),DataFrame(transpose(value.(mdict[m].ext[:expressions][:profit]).data),:auto), delim=";");
        CSV.write(joinpath(home_dir,string("Results_",data["nReprDays"],"_repr_days"),string("Costs",agents[:ps][mm],".csv")),DataFrame(transpose(value.(mdict[m].ext[:expressions][:cost]).data),:auto), delim=";");

        if m == "Edemand"
            cap[mm] = 0
            CVAR[mm] = value.(mdict[m].ext[:expressions][:CVAR])
            VAR[mm] = value.(mdict[m].ext[:variables][:α])
        else
            cap[mm] = value.(mdict[m].ext[:variables][:cap])
            CVAR[mm] = value.(mdict[m].ext[:expressions][:CVAR])
            VAR[mm] = value.(mdict[m].ext[:variables][:α])
            CSV.write(joinpath(home_dir,string("Results_",data["nReprDays"],"_repr_days"),string("Revenues",agents[:ps][mm],".csv")),DataFrame(transpose(value.(mdict[m].ext[:expressions][:revenue]).data),:auto), delim=";");
        end
        for jy in 1:data["nYears"]
            CSV.write(joinpath(home_dir,string("Results_",data["nReprDays"],"_repr_days"),string("Generation",agents[:ps][mm],"_",jy,".csv")), DataFrame(results["g"][agents[:ps][mm]][end][:,:,jy],:auto), delim=";");
        end

        mm = mm + 1
    end
    CSV.write(joinpath(home_dir,string("Results_",data["nReprDays"],"_repr_days"),string("capacity_ps.csv")), DataFrame(transpose(cap),:auto), delim=";",header=string.("CAP_",agents[:ps]));
    CSV.write(joinpath(home_dir,string("Results_",data["nReprDays"],"_repr_days"),string("CVAR_ps.csv")), DataFrame(transpose(CVAR),:auto), delim=";",header=string.("CVAR_",agents[:ps]));
    CSV.write(joinpath(home_dir,string("Results_",data["nReprDays"],"_repr_days"),string("VAR_ps.csv")), DataFrame(transpose(VAR),:auto), delim=";",header=string.("VAR_",agents[:ps]));

    ENS = EOM["D"] + results["g"]["Edemand"][end]
    ENS_H2 = H2["D"] + results["gH"]["H2demand"][end]
    for jy in 1:data["nYears"]
        CSV.write(joinpath(home_dir,string("Results_",data["nReprDays"],"_repr_days"),string("ENS_el_",jy,".csv")), DataFrame(ENS[:,:,jy],:auto), delim=";");


        # Elastic demand
        CSV.write(joinpath(home_dir,string("Results_",data["nReprDays"],"_repr_days"),string("Elastic_demand_el_",jy,".csv")), DataFrame(results["g_ela"][end][:,:,jy],:auto), delim=";");
        CSV.write(joinpath(home_dir,string("Results_",data["nReprDays"],"_repr_days"),string("Elastic_demand_H2_",jy,".csv")), DataFrame(results["gH2_ela"][end][:,:,jy],:auto), delim=";");

        # Electricity prices
        CSV.write(joinpath(home_dir,string("Results_",data["nReprDays"],"_repr_days"),string("electricity_price_",jy,".csv")), DataFrame(results["λ"]["EOM"][end][:,:,jy],:auto), delim=";");

        # Hydrogen prices
        CSV.write(joinpath(home_dir,string("Results_",data["nReprDays"],"_repr_days"),string("hydrogen_price_",jy,".csv")), DataFrame(results["λ"]["H2"][end][:,:,jy],:auto), delim=";");

        # Hydrogen sector
        CSV.write(joinpath(home_dir,string("Results_",data["nReprDays"],"_repr_days"),string("SOC",agents[:h2s][3],"_",jy,".csv")), DataFrame(results["SOC"][end][:,:,jy],:auto), delim=";");
        CSV.write(joinpath(home_dir,string("Results_",data["nReprDays"],"_repr_days"),string("SOC_AD_0",agents[:h2s][3],"_",jy,".csv")), DataFrame(results["SOC_AD_0"][end][:,jy]',:auto), delim=";");
        CSV.write(joinpath(home_dir,string("Results_",data["nReprDays"],"_repr_days"),string("Charge",agents[:h2s][3],"_",jy,".csv")), DataFrame(results["chH"][end][:,:,jy],:auto), delim=";");
        CSV.write(joinpath(home_dir,string("Results_",data["nReprDays"],"_repr_days"),string("Discharge",agents[:h2s][3],"_",jy,".csv")), DataFrame(results["dhH"][end][:,:,jy],:auto), delim=";");
        CSV.write(joinpath(home_dir,string("Results_",data["nReprDays"],"_repr_days"),string("demand_H2s_",jy,".csv")), DataFrame(H2["D"][:,:,jy],:auto), delim=";");
    end

    # Hydrogen sector
    mm = 1
    capH = zeros(length(agents[:h2s]))
    CVAR_H = zeros(length(agents[:h2s]))
    VAR_H = zeros(length(agents[:h2s]))
    volH = zeros(1) # set at one because there is only 1 storage asset, can be increased
    
    for m in agents[:h2s]
 
        CSV.write(joinpath(home_dir,string("Results_",data["nReprDays"],"_repr_days"),string("Diff_VAR",agents[:h2s][mm],".csv")),DataFrame(transpose(value.(mdict[m].ext[:variables][:u]).data),:auto), delim=";");
        CSV.write(joinpath(home_dir,string("Results_",data["nReprDays"],"_repr_days"),string("Profit",agents[:h2s][mm],".csv")),DataFrame(transpose(value.(mdict[m].ext[:expressions][:profit]).data),:auto), delim=";");
        CSV.write(joinpath(home_dir,string("Results_",data["nReprDays"],"_repr_days"),string("Costs",agents[:h2s][mm],".csv")),DataFrame(transpose(value.(mdict[m].ext[:expressions][:cost]).data),:auto), delim=";");

        if m == "H2demand"
            capH[mm] = 0
            CVAR_H[mm] = value.(mdict[m].ext[:expressions][:CVAR])
            VAR_H[mm] = value.(mdict[m].ext[:variables][:α])
        else
            capH[mm] = value.(mdict[m].ext[:variables][:capH])
            CVAR_H[mm] = value.(mdict[m].ext[:expressions][:CVAR])
            VAR_H[mm] = value.(mdict[m].ext[:variables][:α])
            
            if m == "H2storage"
                volH[1] = value.(mdict[m].ext[:variables][:volH])
            end
        end
        for jy in 1:data["nYears"]
            CSV.write(joinpath(home_dir,string("Results_",data["nReprDays"],"_repr_days"),string("H2Generation",agents[:h2s][mm],"_",jy,".csv")), DataFrame(results["h2"][agents[:h2s][mm]][end][:,:,jy],:auto), delim=";");
        end

        mm = mm + 1  
    end
    CSV.write(joinpath(home_dir,string("Results_",data["nReprDays"],"_repr_days"),string("capacity_h2s.csv")), DataFrame(transpose(capH),:auto), delim=";",header=string.("CAP_",agents[:h2s]));
    CSV.write(joinpath(home_dir,string("Results_",data["nReprDays"],"_repr_days"),string("CVAR_h2s.csv")), DataFrame(transpose(CVAR_H),:auto), delim=";",header=string.("CVAR_",agents[:h2s]));
    CSV.write(joinpath(home_dir,string("Results_",data["nReprDays"],"_repr_days"),string("VAR_h2s.csv")), DataFrame(transpose(VAR_H),:auto), delim=";",header=string.("VAR_",agents[:h2s]));
    CSV.write(joinpath(home_dir,string("Results_",data["nReprDays"],"_repr_days"),string("storage_volume.csv")), DataFrame(transpose(volH),:auto), delim=";");

    # Capacity market
    CSV.write(joinpath(home_dir,string("Results_",data["nReprDays"],"_repr_days"),string("capacity_price.csv")), DataFrame(transpose([results["λ"]["CM"][end]]),:auto), delim=";");

    mm = 1
    cap_cm = zeros(length(agents[:cm]))
    for m in agents[:cm]
        cap_cm[mm] = value.(mdict[m].ext[:variables][:cap_cm])
        mm = mm + 1
    end
    CSV.write(joinpath(home_dir,string("Results_",data["nReprDays"],"_repr_days"),string("capacity_offered_cm.csv")), DataFrame(transpose(cap_cm),:auto), header=string.("CAP_",agents[:cm]));

    mm = 1
    capH_cm = zeros(length(agents[:hcm]))
    for m in agents[:hcm]
        capH_cm[mm] = value.(mdict[m].ext[:variables][:capH_cm])
        mm = mm + 1
    end
    CSV.write(joinpath(home_dir, string("Results_", data["nReprDays"], "_repr_days"), string("capacity_offered_H2cm.csv")), DataFrame(transpose(capH_cm), :auto), header=string.("CAP_H_", agents[:hcm]))

    # Net profits

    obj = Dict()
    for m in agents[:all]
        obj[m] = JuMP.value(objective_value(mdict[m]))
    end
    CSV.write(joinpath(home_dir,string("Results_",data["nReprDays"],"_repr_days"),string("objective.csv")), DataFrame(agent = collect(keys(obj)), weigthed_obj = collect(values(obj))), delim=";");

    # Chck convergence
    #CSV.write(joinpath(home_dir, string("Results_", data["nReprDays"], "_repr_days"), string("PrimalResidualEOM.csv")), DataFrame(ADMM["Residuals"]["Primal"]["EOM"][:]', :auto), delim=";")
    #CSV.write(joinpath(home_dir, string("Results_", data["nReprDays"], "_repr_days"), string("PrimalResidualH2.csv")), DataFrame(ADMM["Residuals"]["Primal"]["H2"][:]', :auto), delim=";")
    #CSV.write(joinpath(home_dir, string("Results_", data["nReprDays"], "_repr_days"), string("DualResidualEOM.csv")), DataFrame(ADMM["Residuals"]["Dual"]["EOM"][:]', :auto), delim=";")
    #CSV.write(joinpath(home_dir, string("Results_", data["nReprDays"], "_repr_days"), string("DualResidualH2.csv")), DataFrame(ADMM["Residuals"]["Dual"]["H2"][:]', :auto), delim=";")

    plot(1:150, [ADMM["Residuals"]["Primal"]["EOM"][:] ADMM["Residuals"]["Primal"]["H2"][:] ADMM["Residuals"]["Dual"]["EOM"][:] ADMM["Residuals"]["Dual"]["H2"][:]])

end
