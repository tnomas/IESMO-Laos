# -*- coding: utf-8 -*-

from pyomo.environ import ConcreteModel, Param, Var, Constraint, Objective
from pyomo.environ import minimize, NonNegativeReals
from pyomo.opt import SolverFactory
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from datetime import datetime

# =============================================================================
# --------Initialize-----------
# =============================================================================

start = datetime.now()
data = pd.read_csv('data.csv', header=0)
parameter = pd.read_csv('parameter.csv', index_col='name', header=0)
m = ConcreteModel()
m.T = [i for i in range(0, 8760)]
print('Wait...')

# =============================================================================
# ------------Variables------------
# =============================================================================

# Incentives
Avoid_EmptyStorage = -1
Avoid_Valve = 99999
Avoid_RE = 99999999  # incentive to Avoid_RE

# Wind
C_Wind = parameter.loc['C_Wind', 'value']  # Cost [Wind €/MW/a]
Factor_Wind = data['windfactor']  # Wind factor (t) []
Lt_Wind = parameter.loc['Lt_Wind', 'value']  # Lifetime of wind turbines [a]

# PV
C_Pv = parameter.loc['C_Pv', 'value']  # Cost PV [€/MWp/a]
Factor_Pv = data['pvfactor']  # Radiation factor (t) []
Lt_Pv = parameter.loc['Lt_Pv', 'value']  # Lifetime of PV modules [a]

# Dam
P_Dam = parameter.loc['P_Dam', 'value']  # Max possible turbine power[MW]
Factor_Dam = parameter.loc['Factor_Dam', 'value']  # water usage per MW[m^3/MW]
Sto_Size = parameter.loc['Sto_Size', 'value']  # Storage capacity dam [m^3]
Sto_Start = Sto_Size / 2  # Storage volume @ hour 0  [m^3]
Factor_Inflow = parameter.loc['Factor_inflow', 'value']  # factor to 2017 [%]
WaterInflow = 0.01 * Factor_Inflow * \
              data['riverflow_first'].interpolate(method='linear')
              #  Water Inflow River (t) [m^3]

# Demand
Growth_Demand_E = parameter.loc['Growth_Demand_E', 'value']  # Growth in 2050
Factor_Demand_E = data['energy_d_normed']  # Energy demand factor (t) []
Total_Demand_E = parameter.loc['Total_Demand_E', 'value']  # Total [MWh/a]
Demand_E = [(Growth_Demand_E * Total_Demand_E *
             Factor_Demand_E[i]) for i in m.T]

Growth_Demand_W = parameter.loc['Growth_Demand_W', 'value']  # Growth in 2050
Total_Demand_W = parameter.loc['Total_Demand_W', 'value']  # Total [m^3/a]
Factor_Demand_W = data['water_d_normed']  # Water Demand (t) [m^3]
Demand_W = [(Growth_Demand_W*Total_Demand_W * Factor_Demand_W[i]) for i in m.T]

# =============================================================================
# ------------MODEL CONSTRUCTION------------
# =============================================================================

# Parameter
m.C_Wind = Param(initialize=C_Wind)
m.C_Pv = Param(initialize=C_Pv)
m.Sto_Size = Param(initialize=Sto_Size)
m.Sto_Start = Param(initialize=Sto_Start)
m.Avoid_Valve = Param(initialize=Avoid_Valve)
m.Avoid_RE = Param(initialize=Avoid_RE)
m.Avoid_EmptyStorage = Param(initialize=Avoid_EmptyStorage)

# Variablen
m.P_Wind = Var(domain=NonNegativeReals)  # installed MW Wind
m.P_Wind_Usage = Var(m.T, domain=NonNegativeReals)  # used Wind @ each h
m.P_Wind_Over = Var(m.T, domain=NonNegativeReals)  # overproduced Wind @ each h
m.P_Pv = Var(domain=NonNegativeReals)  # installed MW Pv
m.P_Pv_Usage = Var(m.T, domain=NonNegativeReals)  # used PV @ each h
m.P_Pv_Over = Var(m.T, domain=NonNegativeReals)  # overproduced PV @ each h
m.P_Dam = Var(m.T, bounds=(0, P_Dam))  # used Dam @ each h
m.V_Water = Var(m.T, domain=NonNegativeReals)  # used Water for Power @ each h
m.V_Storage = Var(m.T, bounds=(0, Sto_Size))  # Storage lvl @ each h
m.V_Valve = Var(m.T, domain=NonNegativeReals)  # Water release @ each h

# =============================================================================
# ----------Constraint & Objective & Solver------------
# =============================================================================

# ------Objective Function------
def obj_rule(m):
        return(m.C_Wind * m.P_Wind + m.C_Pv * m.P_Pv +
               sum((m.P_Pv_Over[i] + m.P_Wind_Over[i]) for i in m.T) *
               m.Avoid_RE + sum((m.V_Valve[i]) for i in m.T) * m.Avoid_Valve +
               sum((m.V_Storage[i]) for i in m.T) * m.Avoid_EmptyStorage)


m.min = Objective(sense=minimize, rule=obj_rule)


# ------CONSTRAINTS------
# Fulfil Energy Demand
def demand_energy_rule(m, i):
    return (m.P_Wind_Usage[i] + m.P_Pv_Usage[i] + m.P_Dam[i] == Demand_E[i])


m.EnergyDemand = Constraint(m.T, rule=demand_energy_rule)


# Turbine Capacity Limit
def turbine_limitation_rule(m, i):
    return(m.V_Water[i] == m.P_Dam[i] * Factor_Dam)


m.TurbineLmt = Constraint(m.T, rule=turbine_limitation_rule)


#  Split installed Wind into Used and Overproduced Wind
def wind_rule(m, i):
    return(m.P_Wind_Usage[i] + m.P_Wind_Over[i] == m.P_Wind * Factor_Wind[i])


m.WindP = Constraint(m.T, rule=wind_rule)


#  Split installed PV into Used and Overproduced PV
def pv_rule(m, i):
    return(m.P_Pv_Usage[i] + m.P_Pv_Over[i] == m.P_Pv * Factor_Pv[i])


m.MaxPvP = Constraint(m.T, rule=pv_rule)


# Storage Calculation
def storage_rule(m, i):
    if i == 0:
        return(m.V_Storage[i] == m.Sto_Start - WaterInflow[i] -
               m.V_Water[i] - Demand_W[i] - m.V_Valve[i])
    else:
        return(m.V_Storage[i] == m.V_Storage[i-1] + WaterInflow[i] -
               m.V_Water[i] - Demand_W[i] - m.V_Valve[i])


m.StorageCalc = Constraint(m.T, rule=storage_rule)


# Storage @ hour 8759 = Storage @ hour 0
def dam_end_rule(m, i):
    return(m.V_Storage[8759] == m.Sto_Start)


m.DamEndLvl = Constraint(m.T, rule=dam_end_rule)


# ------SOLVER------
solver = ('gurobi')
opt = SolverFactory(solver)
m.write('optimization_problem.lp', io_options={'symbolic_solver_labels': True})
results = opt.solve(m, tee=True)

# =============================================================================
# ------------Results------------
# =============================================================================

P_Wind_Inst = [m.P_Wind.value for i in m.T]
P_Wind_Out = [(Factor_Wind[i]*m.P_Wind.value) for i in m.T]
P_Wind_Usage = [m.P_Wind_Usage[i].value for i in m.T]
P_Wind_Over = [m.P_Wind_Over[i].value for i in m.T]
P_Pv_Inst = [m.P_Pv.value for i in m.T]
P_Pv_Out = [(Factor_Pv[i]*m.P_Pv.value) for i in m.T]
P_Pv_Usage = [m.P_Pv_Usage[i].value for i in m.T]
P_Pv_Over = [m.P_Pv_Over[i].value for i in m.T]
P_Dam_Out = [(m.V_Water[i].value/Factor_Dam) for i in m.T]
Used_V_Water = [m.V_Water[i].value for i in m.T]
Storage = [m.V_Storage[i].value for i in m.T]
WaterBalance = [(WaterInflow[i] - m.V_Water[i].value -
                 Demand_W[i]) for i in m.T]
Valve = [m.V_Valve[i].value for i in m.T]

Results = pd.DataFrame({"Hour": pd.Series(m.T),
                        "Energy Demand": pd.Series(Demand_E),
                        "Water Demand": pd.Series(Demand_W),
                        "Factor Wind": pd.Series(Factor_Wind),
                        "Wind Output": pd.Series(P_Wind_Out),
                        "Wind Usage": pd.Series(P_Wind_Usage),
                        "Wind Over": pd.Series(P_Wind_Over),
                        "Factor PV": pd.Series(Factor_Pv),
                        "PV Over": pd.Series(P_Pv_Over),
                        "PV Usage": pd.Series(P_Pv_Usage),
                        "PV Output": pd.Series(P_Pv_Out),
                        "Dam Output": pd.Series(P_Dam_Out),
                        "Turbine Water": pd.Series(Used_V_Water),
                        "Storage": pd.Series(Storage),
                        "River": pd.Series(WaterInflow),
                        "Water Balance": pd.Series(WaterBalance),
                        "Valve": pd.Series(Valve),
                        "Installed Wind": pd.Series(P_Wind_Inst),
                        "Installed PV": pd.Series(P_Pv_Inst)},
                       columns=['Hour', 'Energy Demand', 'Water Demand',
                                'Factor Wind', 'Wind Output', 'Wind Usage',
                                'Wind Over', 'Factor PV', 'PV Output',
                                'PV Usage', 'PV Over', 'Dam Output',
                                'Storage', 'Turbine Water', 'River',
                                'Water Balance', 'Valve', 'Installed Wind',
                                'Installed PV'])

Results = (np.round(Results, decimals=2))
Results.to_csv('output.csv', sep=";", decimal=",")

# =============================================================================
# ------------Plotting------------
# =============================================================================

Graph_Start = 0
Graph_End = 8760

ResultsGraph = Results[Graph_Start:Graph_End]
demand = ResultsGraph['Energy Demand']
y1 = ResultsGraph['Dam Output']
y2 = y1 + ResultsGraph['PV Usage']
y3 = y2 + ResultsGraph['Wind Usage']
y4 = ResultsGraph['Storage']
y5 = ResultsGraph['Valve']

fig = plt.figure()
ax = fig.add_subplot(211)
ax.fill_between(ResultsGraph.Hour, 0, y1, label='Dam', color='g', alpha=.7)
ax.fill_between(ResultsGraph.Hour, y1, y2, label='PV', color='y', alpha=.7)
ax.fill_between(ResultsGraph.Hour, y2, y3, label='Wind', color='b', alpha=.7)
ax.plot(ResultsGraph.Hour, demand, label='Energy Demand', color='black')
ax.yaxis.set_label_position('left')
ax.set_ylabel('[MW]')
plt.legend(fontsize='small', loc='upper center', bbox_to_anchor=(0.6, 1.3),
           fancybox=True, ncol=4)

ax2 = fig.add_subplot(212)
ax2.fill_between(ResultsGraph.Hour, 0, y4, label='Storage',
                 color='orange', alpha=.7)
ax2.plot(ResultsGraph.Hour, y5, label='Valve', color='blue', alpha=.7)
ax2.set_ylabel('[m^3]')
plt.xlabel('Hours')
plt.legend(fontsize='small', loc='upper center', bbox_to_anchor=(0.81, -0.2),
           fancybox=True, ncol=4)

plt.savefig('graphical_output/energy.pdf', dpi=150)

ResultsGraph.to_csv('selection.csv', sep=";", decimal=",")
ResultsGraph.to_csv('alternative.csv')

# =============================================================================
# ------------OUTPUT------------
# =============================================================================

print(np.round(Results, decimals=2))

duration = datetime.now() - start
print("Problem solved in", str(duration), "using", solver, "as solver.\n"
      "Results saved in 'output.csv'")

InstalledWind = round(m.P_Wind.value, 2)
InstalledPV = round(m.P_Pv.value, 2)
print("\n----Installed Capacity---- \n"
      "Dam:", P_Dam, "MW \n"
      "Wind:", InstalledWind, "MW \n"
      "PV:", InstalledPV, "MW  \n"
      "\n----COST----\n"
      "Wind:", "{:0,.2f}".format(C_Wind * InstalledWind * Lt_Wind), "€ total",
      "|", C_Wind, "€/MW \n"
      "PV:", "{:0,.2f}".format(C_Pv * InstalledPV*Lt_Pv), "€ total |",
      C_Pv, "€/MW \n"
      "Total:", "{:0,.2f}".format(C_Wind * InstalledWind * Lt_Wind + C_Pv *
                                  InstalledPV * Lt_Pv), "€ \n")
