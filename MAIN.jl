## Topic: Hydrogen market for long-term storage
# Author: Kenneth Bruninx
# Last update: January 2023

## 0. Set-up code
println(string("#################                    Start simulation                         ###################"))
println("        ")

# HPC or not?
HPC = "NA" # NA, DelftBlue or ThinKing

# Home directory
const home_dir = @__DIR__

if HPC == "DelftBlue"  # only for running this on DelftBlue -- not relevant for now
    ENV["GRB_LICENSE_FILE"] = "/home/aberdin/gurobi.lic"
    ENV["GUROBI_HOME"] = "/scratch/aberdin/gurobi952/linux64"
    println(string("Number of threads: ", Threads.nthreads()))
end

# Install
#import Pkg
#Pkg.add("Gurobi")
#Pkg.status()
#Pkg.add("DataFrames")
#Pkg.add("CSV")
#Pkg.add("YAML")
#Pkg.add("DataStructures")
#Pkg.add("ProgressBars")
#Pkg.add("Printf")
#Pkg.add("TimerOutputs")
#Pkg.add("ArgParse")

# Include packages 
using JuMP # Optimization packages
using Gurobi
using DataFrames, CSV, YAML, DataStructures # dataprocessing
using ProgressBars, Printf # progress bar
using TimerOutputs # profiling 
using Base.Threads: @spawn
using Base: split
using ArgParse # Parsing arguments from the command line

# Gurobi environment to suppress output
println("Define Gurobi environment...")
println("        ")
const GUROBI_ENV = Gurobi.Env()
# set parameters:
GRBsetparam(GUROBI_ENV, "OutputFlag", "0")
println("        ")

# Include functions
include(joinpath(home_dir, "Source", "define_common_parameters.jl"))
include(joinpath(home_dir, "Source", "define_H2S_parameters.jl"))
include(joinpath(home_dir, "Source", "define_ps_parameters.jl"))
include(joinpath(home_dir, "Source", "define_EOM_parameters.jl"))
include(joinpath(home_dir, "Source", "define_H2_parameters.jl"))
include(joinpath(home_dir, "Source", "build_ps_agent.jl"))
include(joinpath(home_dir, "Source", "build_H2S_agent.jl"))
include(joinpath(home_dir, "Source", "define_results.jl"))
include(joinpath(home_dir, "Source", "ADMM.jl"))
include(joinpath(home_dir, "Source", "ADMM_subroutine.jl"))
include(joinpath(home_dir, "Source", "solve_ps_agent.jl"))
include(joinpath(home_dir, "Source", "solve_H2S_agent.jl"))
include(joinpath(home_dir, "Source", "update_rho.jl"))
include(joinpath(home_dir, "Source", "save_results.jl"))

# Data 
data = YAML.load_file(joinpath(home_dir, "Input", "overview_data.yaml"))

ts = Dict()
order_matrix = Dict()
repr_days = Dict()
# years = Dict(1 => 2021) # deterministic
# years = Dict(1 => 2017, 2 => 2018, 3 => 2019, 4 => 2020, 5 => 2021, 6 => 2022) # stochastic
# years = Dict(1 => 2021, 2 => 20211) # validation stochastic
# years = Dict(1 => 2017, 2 => 2018, 3 => 2019, 4 => 2020) # validation CVAR
years = Dict(1 => "2017", 2 => "2017_Hhigh", 3 => "2017_Hlow", 4 => "2018", 5 => "2018_Hhigh", 6 => "2018_Hlow", 7 => "2019", 8 => "2019_Hhigh", 9 => "2019_Hlow", 
10 => "2020", 11 => "2020_Hhigh", 12 => "2020_Hlow", 13 => "2021", 14 => "2021_Hhigh", 15 => "2021_Hlow", 16 => "2022", 17 => "2022_Hhigh", 18 => "2022_Hlow" )

for yr in keys(years)
    ts[yr] = CSV.read(joinpath(home_dir, "Input", "timeseries", string("timeseries_", years[yr], ".csv")), delim=",", DataFrame)
    order_matrix[yr] = CSV.read(joinpath(home_dir, "Input", string("output_",years[yr]), "ordering_variable.csv"), delim=",", DataFrame)
    #order_matrix[yr] = CSV.read(joinpath(home_dir, "Input", "output_2021", "ordering_variable.csv"), delim=",", DataFrame)
    repr_days[yr] = CSV.read(joinpath(home_dir, "Input", string("output_",years[yr]), "decision_variables_short.csv"), delim=",", DataFrame)
    #repr_days[yr] = CSV.read(joinpath(home_dir, "Input", "output_2021", "decision_variables_short.csv"), delim=",", DataFrame)
end

# Create folder for results
if isdir(joinpath(home_dir, string("Results_", data["General"]["nReprDays"], "_repr_days"))) != 1
    mkdir(joinpath(home_dir, string("Results_", data["General"]["nReprDays"], "_repr_days")))
end

println("Including all required input data: done")
println("   ")

## 2. Initiate models for representative agents 
agents = Dict()
agents[:ps] = [id for id in keys(data["PowerSector"])]
agents[:h2s] = [id for id in keys(data["HydrogenSector"])]
agents[:all] = union(agents[:ps], agents[:h2s])
# Different markets - to be completed based on the agents
agents[:eom] = []
agents[:h2] = []
agents[:cm] = []
agents[:hcm] = []
mdict = Dict(i => Model(optimizer_with_attributes(() -> Gurobi.Optimizer(GUROBI_ENV))) for i in agents[:all])

## 3. Define parameters for markets and representative agents

EOM = Dict{String,Any}()
H2 = Dict{String,Any}()
CM = Dict{String,Any}()
HCM = Dict{String,Any}()

# Parameters/variables EOM
define_EOM_parameters!(EOM, merge(data["General"], data["EOM"]), ts, repr_days)

# Parameters/variables Hydrogen Market
define_H2_parameters!(H2, merge(data["General"], data["H2"]), ts, repr_days)

# Parameters/variables Capacity Market
CM["D"] =  12.525  # GW
HCM["D"] = 0 # GW

# Parameters agents
for m in agents[:ps]
    define_common_parameters!(m, mdict[m], merge(data["General"], data["ADMM"], data["PowerSector"][m]), ts, repr_days, agents)     # Parameters common to all agents
    define_ps_parameters!(m, mdict[m], merge(data["General"], data["PowerSector"][m]), ts, repr_days)                               # Power sector
end
for m in agents[:h2s]
    define_common_parameters!(m, mdict[m], merge(data["General"], data["ADMM"], data["HydrogenSector"][m]), ts, repr_days, agents)  # Parameters common to all agents 
    define_H2S_parameters!(m, mdict[m], merge(data["General"], data["HydrogenSector"][m]), repr_days, order_matrix)                 # Hydrogen sector
end

# Calculate number of agents in each market
EOM["nAgents"] = length(agents[:eom])
H2["nAgents"] = length(agents[:h2])
CM["nAgents"] = length(agents[:cm])
HCM["nAgents"] = length(agents[:hcm])

println("Inititate model, sets and parameters: done")
println("   ")

## 4. Build models
for m in agents[:ps]
    build_ps_agent!(m, mdict[m], EOM)
end
for m in agents[:h2s]
    build_h2s_agent!(m, mdict[m], H2)
end

println("Build model: done")
println("   ")

## 5. ADMM proces to calculate equilibrium
println("Find equilibrium solution...")
println("   ")
println("(Progress indicators on primal residuals, relative to tolerance: <1 indicates convergence)")
println("   ")

results = Dict()
ADMM = Dict()
TO = TimerOutput()
define_results!(merge(data["General"], data["ADMM"]), results, ADMM, agents, EOM, H2)             # initialize structure of results, only those that will be stored in each iteration
ADMM!(results, ADMM, EOM, H2, mdict, agents, data, TO)                                             # calculate equilibrium 
ADMM["walltime"] = TimerOutputs.tottime(TO) * 10^-9 / 60                                       # wall time 

println(string("Done!"))
println(string("        "))
println(string("Required iterations: ", ADMM["n_iter"]))
println(string("Required walltime: ", ADMM["walltime"], " minutes"))
println(string("        "))

## 6. Postprocessing and save results 
save_results(mdict, EOM, H2, ADMM, results, merge(data["General"], data["ADMM"], data["H2"]), agents)
YAML.write_file(joinpath(home_dir, string("Results_", data["General"]["nReprDays"], "_repr_days"), string("TimerOutput.yaml")), TO)

println("Postprocessing & save results: done")
println("   ")

println(string("##############################################################################################"))