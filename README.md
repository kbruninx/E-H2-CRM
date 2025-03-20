This code was developed to study the interaction between the power sector and the hydrogen sector in the context of interdependents energy-only electricity market, hydrogen market and potential capacity remuneration mechanisms for both the energy vectors. It calculates an equilibrium between a set of representative price-taking agents on these markets. It allows to enforce hydrogen and renewable capacity targets. It employs an iterative price-search algorithm based on ADMM to calculate this equilibrium iteratively.

The research code accompanies the following paper: "Risk-aversion and Capacity Remuneration Mechanisms in Electricity-Hydrogen Markets". 

## Installation, hardware & software requirements
### Installation
After downloading the repository, no specific installation is required. The user should, however, have the following software installed to be able to execute the program:
- Julia (https://julialang.org/downloads/)
- Gurobi (https://www.gurobi.com/) and have a license for this solver. If the user doesn't have access to Gurobi, any alternative (open source) solver capable of solving quadratic programming problems can be used. 

The program can be executed in the Julia terminal. Any text editor can be used to modify the code. However, it may be more convenient to work with a tool such as Visual Studio Code (https://code.visualstudio.com/).

If the user does not have any of these programs installed, the installation should take approximately one hour.

### Software requirements: 
This code has been developed using Julia v1.8. The solver used is Gurobi v.9.5.2.

The following Julia packages are required:
- JuMP
- Gurobi
- DataFrames
- CSV
- YAML
- DataStructures
- ProgressBars
- Printf
- TimerOutputs
- ArgParse

If the user does not have these programs installed, the installation should take less than one hour.

### Hardware requirements 
No specific hardware is required. Depending on the configuration (number of agents and markets considered), computational effort may significantly increase.

## Running the program
### Input
The file "overview_data.yaml" contains a number of input parameters that are common to all scenarios. These includes general parameters for the simulation, for the algorithm, and market and technologies parameters.

In particular, to obtain the results presented in the paper the user can modify:

1. The weight of risk measure "gamma" at line 8, in order to choose between the risk-averse (gamma = 0.5) and risk-neutral (gamma = 1) setups.
2. The risk aversion parametrization "beta" at line 9, in order to set the degree of risk aversion considered by the CVAR measure (0.2 <= beta <= 1).
3. The switch "sigmaCM" for the electricity capacity market at line 10, in order to choose if a capacity market for electricity is active (sigmaCM = 1) or not (sigmaCM = 0).
4. The switch "sigmaHCM" for the hydrogen capacity market at line 11, in order to choose if a capacity market for hydrogen is active (sigmaHCM = 1) or not (sigmaHCM = 0).

The ADMM parameters may require tuning depending on the degree of risk aversion chosen.

### Executing the code
The code can be run by executing the "MAIN.jl" file.

### Output & Postprocessing
Running the code will generate the following output files for each run, located in the "Results_8_repr_days":

    1.   "capacity_h2s.csv": hydrogen installed capacity for each technology (GW)
    2.   "capacity_offered_cm.csv": annual power capacity offered in the capacity market for each participating technology (GW)
    3.   "capacity_offered_H2cm.csv": annual hydrogen capacity offered in the capacity market for each participating technology (GW)
    4.   "capacity_price.csv": electricity capacity market price (M€/GW)
    5.   "capacity_ps.csv": power installed capacity for each technology (GW)
    6.   "ChargeH2storage_X.csv" (X = scenario number): hourly hydrogen storage charge for each representative day (MWh_H2)
    7.   "CostsY.csv" (Y = agent): sum of costs of agent Y for each scenario (M€)
    8.   "CVAR_h2s.csv": value of CVAR for each hydrogen agent (M€)
    9.   "CVAR_ps.csv": value of CVAR for each power agent (M€)
    10.  "demand_H2s_X.csv" (X = scenario number): hourly hydrogen demand for each representative day (MWh_H2)
    11.  "demand_ps_X.csv" (X = scenario number): hourly power demand for each representative day (MWh)
    12.  "Diff_VARY.csv" (Y = agent): difference in profit with the VAR for each scenario (M€)
    13.  "DischargeH2storage_X.csv" (X = scenario number): hourly hydrogen storage discharge for each representative day (MWh_H2)
    14.  "Elastic_demand_el_X.csv" (X = scenario number): hourly flexible power demand component for each representative day (MWh_H2)
    15.  "Elastic_demand_H2_X.csv" (X = scenario number): hourly flexible hydrogen demand component for each representative day (MWh_H2)
    16.  "electricity_price_X.csv" (X = scenario number): hourly electricity price for each representative day (€/MWh)
    17.  "ENS_el_X.csv" (X = scenario number): hourly electricity not served with respect to nominal electricity demand for each representative day (MWh)
    18.  "ENS_H2_X.csv" (X = scenario number): hourly  hydrogen not served with respect to nominal hydrogen demand for each representative day (MWh_H2)
    19.  "GenerationY_X" (X = scenario number, Y = agent): hourly electricity generation (or demand, when negative) of the Y agent for each representative day (MWh)
    20.  "H2GenerationY_X" (X = scenario number, Y = agent): hourly hydrogen generation (or demand, when negative) of the Y agent for each representative day (MWh_H2)
    21.  "hydrogen_price_X.csv" (X = scenario number): hourly hydrogen price for each representative day (€/MWh_H2)
    22.  "objective.csv": value of the objective function for each agent (M€)
    23.  "ProfitY.csv" (Y = agent): expected profit of agent Y for each scenario (M€)
    24.  "RevenuesY.csv" (Y = agent): sum of revenues of agent Y for each scenario (M€)
    25.  "SOC_AD_0H2storage_X.csv" (X = scenario number): initial state of charge of hydrogen storage for each day (TWh)
    26.  "SOCH2storage_X.csv" (X = scenario number): hourly state of charge of hydrogen storage for each representative day (TWh)
    27.  "storage_volume.csv": hydrogen storage energy capacity (TWh)
    28.  "VAR_h2s.csv": value-at-risk of each hydrogen market agent (M€)
    29.  "VAR_ps.csv": value-at-risk of each electricity market agent (M€)

## License
The software is made available under the MIT license (https://opensource.org/licenses/MIT).
 
## Contributors
A. Berdin (alessio.berdin1999@gmail.com) and K. Bruninx (k.bruninx@tudelft.nl)
