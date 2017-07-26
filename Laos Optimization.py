#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul 13 17:52:14 2017
@author: nilshoffmann
"""

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
m = ConcreteModel()
m.T = [i for i in range(0, 8760)]
print('Wait...')

# =============================================================================
# ------------Variables------------
# =============================================================================

# Incentives
Prio_Wind = 3  # incentive to use P_Dam over P_Wind
Prio_PV = 3  # incentive to use P_Dam over P_Pv
Prio_Dam = 1  # incentive to use P_Dam first
Prio_Sto = 0.000000000000001
Avoid = 99999999  # incentive to avoid

# Wind
C_Wind = 77350  # Cost [Wind €/MW/a]
Factor_Wind = data['windfactor']  # Wind factor (t) []
Lt_Wind = 20  # Lifetime of wind turbines [a]

# PV
C_Pv = 37120  # Cost PV [€/MWp/a]
Factor_Pv = data['pvfactor']  # Radiation factor (t) []
Lt_Pv = 22  # Lifetime of PV modules [a]

# Dam
P_Dam = 260  # Max possible turbine power[MW]
Factor_Dam = 20547  # Used water/MW turbine power [m^3/MW]
Sto_Size = 6240000  # Storage capacity dam [m^3]
Sto_Start = Sto_Size / 2  # Storage volume @ hour 0  [m^3]
Factor_Inflow = 5  # Waterflow compared to 2017 [%]
WaterInflow = 0.01 * Factor_Inflow * \
              data['riverflow_first'].interpolate(method='linear')
              # Water Inflow River (t) [m^3]

# Demand
Factor_Demand_E = data['energy_d_normed']  # Energy demand factor (t) []
Total_Demand_E = 97236  # Energy demand total [MWh]
Demand_E = [(Total_Demand_E * Factor_Demand_E[i]) for i in m.T]
Demand_W = data['water_d_absolut']  # Water Demand (t) [m^3]

# =============================================================================
# ------------MODEL CONSTRUCTION------------
# =============================================================================

# Parameter
m.C_Wind = Param(initialize=C_Wind)
m.C_Pv = Param(initialize=C_Pv)
m.Sto_Size = Param(initialize=Sto_Size)
m.Sto_Start = Param(initialize=Sto_Start)
m.Prio_Wind = Param(initialize=Prio_Wind)
m.Prio_PV = Param(initialize=Prio_PV)
m.Prio_Dam = Param(initialize=Prio_Dam)
m.Prio_Sto = Param(initialize=Prio_Sto)
m.Avoid = Param(initialize=Avoid)

# Variablen
m.P_Wind = Var(domain=NonNegativeReals)  # installed MW Wind
m.P_Wind_Usage = Var(m.T, domain=NonNegativeReals)
m.P_Wind_Over = Var(m.T, domain=NonNegativeReals)
m.P_Pv = Var(domain=NonNegativeReals)  # installed MW Pv
m.P_Pv_Usage = Var(m.T, domain=NonNegativeReals)
m.P_Pv_Over = Var(m.T, domain=NonNegativeReals)
m.P_Dam = Var(m.T, bounds=(0, P_Dam))
m.P_Water = Var(m.T, domain=NonNegativeReals)  # used Water for Power
m.Sto_Balance = Var(m.T, bounds=(0, Sto_Size))  # Storage lvl @ each hour
m.Valve = Var(m.T, domain=NonNegativeReals)  # Valve to release stored water

# =============================================================================
# ----------Constraint & Objective & Solver------------
# =============================================================================

# ------Objective Function------
def obj_rule(m):
        return(m.C_Wind * m.P_Wind + m.C_Pv * m.P_Pv +
               sum(m.P_Wind_Usage[i] for i in m.T) * m.Prio_Wind +
               sum(m.P_Pv_Usage[i] for i in m.T) * m.Prio_PV +
               sum(m.P_Water[i] for i in m.T)/Factor_Dam * m.Prio_Dam +
               sum((m.P_Pv_Over[i] + m.P_Wind_Over[i]) for i in m.T) * m.Avoid-
               sum((m.Sto_Balance[i] * m.Prio_Sto) for i in m.T))

m.cost = Objective(sense=minimize, rule=obj_rule)


# ------CONSTRAINTS------
# Fulfil Energy Demand
def demand_energy_rule(m, i):
    return (m.P_Wind_Usage[i] + m.P_Pv_Usage[i] + m.P_Dam[i] == Demand_E[i])
m.EnergyDemand = Constraint(m.T, rule=demand_energy_rule)


# Turbine Capacity Limit
def turbine_limitation_rule(m, i):
    return(m.P_Water[i] == m.P_Dam[i] * Factor_Dam)
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
        return(m.Sto_Balance[i] == m.Sto_Start - m.Valve[i])
    else:
        return(m.Sto_Balance[i] == m.Sto_Balance[i-1] + WaterInflow[i] -
               m.P_Water[i] - Demand_W[i] - m.Valve[i])
m.StorageCalc = Constraint(m.T, rule=storage_rule)


# Storage @ hour 8759 = Storage @ hour 0
def dam_end_rule(m, i):
    return(m.Sto_Balance[8759] == m.Sto_Start)
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
P_Dam_Out = [(m.P_Water[i].value/Factor_Dam) for i in m.T]
Used_P_Water = [m.P_Water[i].value for i in m.T]
Storage = [m.Sto_Balance[i].value for i in m.T]
WaterBalance = [(WaterInflow[i] - m.P_Water[i].value -
                 Demand_W[i]) for i in m.T]
Valve = [m.Valve[i].value for i in m.T]

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
                        "Turbine Water": pd.Series(Used_P_Water),
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

AnzeigeAnfang = 0
AnzeigeEnde = 8760

def zero_to_nan(values):
    return [float('nan') if x==0 else x for x in values]

Valve_nan = zero_to_nan(Valve)

ResultsGraph = Results[AnzeigeAnfang:AnzeigeEnde]
demand = ResultsGraph['Energy Demand']
y1 = ResultsGraph['Dam Output']
y2 = y1 + ResultsGraph['PV Usage']
y3 = y2 + ResultsGraph['Wind Usage']
y4 = ResultsGraph['Storage']
y5 = y4+Valve_nan

fig = plt.figure()
ax = fig.add_subplot(211)
ax.fill_between(ResultsGraph.Hour, 0, y1, label='Dam', color='m', alpha=.7)
ax.fill_between(ResultsGraph.Hour, y1, y2, label='PV', color='g', alpha=.7)
ax.fill_between(ResultsGraph.Hour, y2, y3, label='Wind', color='b', alpha=.7)
ax.plot(ResultsGraph.Hour, demand, label='Energy Demand', color='black')
ax.set_ylabel('[MW]')
ax.yaxis.set_label_position('left')
box = ax.get_position()
ax.set_position([box.x0, box.y0 + box.height * 0.1,
                 box.width, box.height * 0.9])

plt.legend(fontsize='small', loc='upper center', bbox_to_anchor=(0.5, -0.15),
           fancybox=True, ncol=5)

ax2 = fig.add_subplot(212)
ax2.fill_between(ResultsGraph.Hour, 0, y4, label='Storage',
                 color='orange', alpha=.7)
ax2.plot(ResultsGraph.Hour, y5, label='Valve', color='blue', alpha=.7)
ax2.set_ylabel('[m^3]')
plt.xlabel('Hours')

plt.savefig('graphical_output/energy.pdf', dpi=150)

ResultsGraph.to_csv('selection.csv', sep=";", decimal=",")

# ==============================================================================
# ------------OUTPUT------------
# ==============================================================================

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
      "|", C_Wind, "€ per MW installed \n"
      "PV:", "{:0,.2f}".format(C_Pv * InstalledPV*Lt_Pv), "€ total |",
      C_Pv, "€ per MW installed \n"
      "Total:", "{:0,.2f}".format((C_Wind * InstalledWind * Lt_Wind + C_Pv *
                     InstalledPV * Lt_Pv), "€ \n"))
