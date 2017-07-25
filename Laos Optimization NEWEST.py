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
import time

# ==============================================================================
# --------Initialize-----------
# ==============================================================================

WaterFlowFactor = 5
start = datetime.now()
data = pd.read_csv('data.csv', header=0)
inflow = data['riverflow_first']
inflow = inflow.interpolate(method = 'linear')
model = ConcreteModel()
model.T = [i for i in range(0, 8760)]

# ==============================================================================
# ------------Variables------------
# ==============================================================================

# Wind
Cwind = 77350  # €/MW/a
Cwind_var = 0.001
FactorWind = data['windfactor']  # Wind Factor @ hour X
LifeTimeWind = 20

# PV
Cpv = 37120  # €/MWp/a
Cpv_var = 0.001
FactorPv = data['pvfactor']  # Radiation Factor @ hour X
LifeTimePv = 25

# Dam
StorageSize = 6240000  # m^3
Pdam = 260  # MW
FactorDam = 20547  # m^3 Water per MW
WaterInflow = 0.01 * WaterFlowFactor * inflow #data['riverflow_absolut']
Cdam = 0  # €/MWh
DamBalance_Start = StorageSize / 2
Expensive = 9999999999

# Demand
DemandFactor_Energy = data['energy_demand_normed']
DemandTotal_Energy = 71737.392  # MWh Total

# Costcalculation
DemandEnergy = [(DemandTotal_Energy * DemandFactor_Energy[i]) for i in model.T]
DemandWater = data['water_demand_absolut']
DemandWater_Factor = 1  # from 0 to 1

print('Wait...')

# ==============================================================================
# ------------MODEL CONSTRUCTION------------
# ==============================================================================
# Parameter
model.Cwind = Param(initialize=Cwind)  # Price per MW Wind €/MWh
model.Cwind_var = Param(initialize=Cwind_var)
model.Cpv = Param(initialize=Cpv)  # Price per MW PV €/MWh
model.Cpv_var = Param(initialize=Cpv_var)
model.Cdam = Param(initialize=Cdam)  # Price of Water €/MWh
model.StorageSize = Param(initialize=StorageSize)  # m^3
model.DamBalance_Start = Param(initialize=DamBalance_Start)
model.Expensive = Param(initialize=Expensive)

# Variablen
model.Pwind = Var(domain=NonNegativeReals)  # installed MW Wind
model.Pwind_Usage = Var(model.T, domain=NonNegativeReals)
model.Pwind_Over = Var(model.T, domain=NonNegativeReals)
model.Ppv = Var(domain=NonNegativeReals)  # installed MW Wind
model.Ppv_Usage = Var(model.T, domain=NonNegativeReals)
model.Ppv_Over = Var(model.T, domain=NonNegativeReals)
model.PowerGeneratingWater = Var(model.T, domain=NonNegativeReals)
model.DamBalance = Var(model.T, bounds=(0, StorageSize))
model.Pdam = Var(model.T, bounds=(0, Pdam))
model.Overflow = Var(model.T, domain=NonNegativeReals)

# ==============================================================================
# ----------Constraint & Objective & Solver------------
# ==============================================================================

# ------Objective Function-------
# €/MW * MW + €/MW + m^3(gesamt) * MWh/m^3 * €/MWh, Ziel: Minimieren

i = 0
def obj_rule(model):
        return(model.Cwind * model.Pwind + model.Cpv * model.Ppv +
               sum(model.Pwind_Usage[i] for i in model.T) * model.Cwind_var +
               sum(model.Ppv_Usage[i] for i in model.T) * model.Cpv_var +
               sum(model.PowerGeneratingWater[i] for i in model.T)/FactorDam * model.Cdam +
               (model.Overflow[i] + model.Ppv_Over[i] + model.Pwind_Over[i]) * model.Expensive)

model.cost = Objective(sense=minimize, rule=obj_rule)

# ----CONSTRAINTS-------
# Fullfil Demand
def DemandEnergy_rule(model, i):
    return (model.Pwind_Usage[i] + model.Ppv_Usage[i] + model.Pdam[i] == DemandEnergy[i])

model.EnergyDemand = Constraint(model.T, rule=DemandEnergy_rule)


# Limitation through Turbine Capacity
def MaxWaterPower_rule(model, i):
    return(model.PowerGeneratingWater[i] == model.Pdam[i] * FactorDam)

model.MaxWaterPower = Constraint(model.T, rule=MaxWaterPower_rule)

def MaxWindPower_rule(model, i):
    return(model.Pwind_Usage[i] + model.Pwind_Over[i] == model.Pwind * FactorWind[i])

model.MaxWindPower = Constraint(model.T, rule=MaxWindPower_rule)

def MaxPvPower_rule(model, i):
    return(model.Ppv_Usage[i] + model.Ppv_Over[i] == model.Ppv * FactorPv[i])

model.MaxPvPower = Constraint(model.T, rule=MaxPvPower_rule)


# Calculate Storage
def WaterUsage_rule(model, i):
    if i == 0:
        return(model.DamBalance[i] == model.DamBalance_Start)
    if i == 8759:
        return(model.DamBalance[8759] == model.DamBalance_Start)
    else:
        return(model.DamBalance[i] == model.DamBalance[i-1] + WaterInflow[i] -
               model.PowerGeneratingWater[i] - DemandWater[i] *
               DemandWater_Factor - model.Overflow[i])

model.WaterDemand = Constraint(model.T, rule=WaterUsage_rule)

# Dont use more Water than in Storage


def GenerationWaterLimit_rule(model, i):
    if i == 0:
        return(model.PowerGeneratingWater[i] + DemandWater[i] *
               DemandWater_Factor - WaterInflow[i] +
               model.Overflow[i] <= model.DamBalance_Start)
    else:
        return(model.PowerGeneratingWater[i] + DemandWater[i] *
               DemandWater_Factor - WaterInflow[i] +
               model.Overflow[i] <= model.DamBalance[i])

model.GenerationWaterLimit = Constraint(model.T, rule=GenerationWaterLimit_rule)


## Same Storage in the End level like @ hour 1
#def DamEndLvl_rule(model, i):
#    return(model.DamBalance[8759] == model.DamBalance_Start)
#
#model.DamEndLvl = Constraint(model.T, rule=DamEndLvl_rule)

# ------SOLVER---------
solver = ('gurobi')
opt = SolverFactory(solver)
model.write('optimization_problem.lp',
         io_options={'symbolic_solver_labels': True})
results = opt.solve(model, tee=True)

# ==============================================================================
# ------------Results------------
# ==============================================================================
Hour = []
WindOutput = []
EnergyDemand = []; WaterDemand = []
PVOutput = []; DamOutput = []
GenerationWaterUsed = []; Storage = []
WaterBalance = []; River = []
Overflow = []
Ppv_Usage = []; Ppv_Over = []
Pwind_Usage = []; Pwind_Over = []
for i in model.T:
    Hour.append(i)
    WindOutput.append(FactorWind[i]*model.Pwind.value)
    PVOutput.append(FactorPv[i]*model.Ppv.value)
    DamOutput.append(model.PowerGeneratingWater[i].value/FactorDam)
    GenerationWaterUsed.append(model.PowerGeneratingWater[i].value)
    WaterBalance.append(WaterInflow[i] - model.PowerGeneratingWater[i].value - DemandWater[i])
    River.append(WaterInflow[i])
    Storage.append(model.DamBalance[i].value)
    Overflow.append(model.Overflow[i].value)
    Ppv_Usage.append(model.Ppv_Usage[i].value)
    Ppv_Over.append(model.Ppv_Over[i].value)
    Pwind_Usage.append(model.Pwind_Usage[i].value)
    Pwind_Over.append(model.Pwind_Over[i].value)

Results = pd.DataFrame({"Hour": pd.Series(Hour), 
                        "Energy Demand": pd.Series(DemandEnergy),
                        "Water Demand": pd.Series(DemandWater),
                        "Wind Output": pd.Series(WindOutput),
                        "PV Output": pd.Series(PVOutput),
                        "Dam Output": pd.Series(DamOutput),
                        "Turbine Water": pd.Series(GenerationWaterUsed),
                        "Factor Wind": pd.Series(FactorWind),
                        "Factor PV": pd.Series(FactorPv),
                        "Storage": pd.Series(Storage),
                        "River": pd.Series(River),
                        "Water Balance": pd.Series(WaterBalance),
                        "Overflow": pd.Series(Overflow),
                        "PV Usage": pd.Series(Ppv_Usage),
                        "PV Over": pd.Series(Ppv_Over),
                        "Wind Usage": pd.Series(Pwind_Usage),
                        "Wind Over": pd.Series(Pwind_Over)},
                        columns=['Hour','Energy Demand','Wind Output','PV Output',
                                 'Dam Output','Factor Wind','Factor PV','Storage',
                                 'Turbine Water','Water Demand','River',
                                 'Water Balance','Overflow','PV Usage',
                                 'PV Over','Wind Usage','Wind Over'])
                        
Results = (np.round(Results, decimals=2))  
Results.to_csv('output.csv', sep = ";", decimal=",")
#==============================================================================
# ------------Plotting------------
#==============================================================================

AnzeigeAnfang = 0 #int(input("Ab welcher Stunde möchtest du die Werte ansehen: "))
AnzeigeEnde = 8760 #int(input("Bis welcher Stunde möchtest du die Werte ansehen: "))

ResultsGraph = Results[AnzeigeAnfang:AnzeigeEnde]
demand = ResultsGraph['Energy Demand']
y1 = ResultsGraph['PV Usage']
y2 = y1 + ResultsGraph['Dam Output']
y3 = y2 + ResultsGraph['Wind Usage']
y4 = ResultsGraph['Storage']
#y5 = y1 + ResultsGraph['Excess']

fig = plt.figure()
ax = fig.add_subplot(211)
ax.fill_between(ResultsGraph.Hour, 0, y1, label = 'PV', color = 'm', alpha=.7)
ax.fill_between(ResultsGraph.Hour, y1, y2, label = 'Dam', color = 'g', alpha=.7)
ax.fill_between(ResultsGraph.Hour, y2, y3, label = 'Wind', color = 'b', alpha=.7)
ax.plot(ResultsGraph.Hour, demand, label = 'Energy Demand', color = 'black')
ax.set_ylabel('[MW]')
ax.yaxis.set_label_position('left')
box = ax.get_position()
ax.set_position([box.x0, box.y0 + box.height * 0.1, box.width, box.height * 0.9])

plt.legend(fontsize = 'small', loc='upper center', bbox_to_anchor=(0.5, -0.15), fancybox=True, ncol=5)

ax2 = fig.add_subplot(212)
ax2.fill_between(ResultsGraph.Hour, y1, y4, label = 'Storage', color = 'orange', alpha=.7)
ax2.set_ylabel('[m^3]')
plt.xlabel('Hours')

plt.savefig('graphical_output/energy.pdf', dpi=150)

ResultsGraph['Overflow'].plot(label='Overflow')
ResultsGraph.to_csv('selection.csv', sep = ";", decimal=",")

#==============================================================================
# ------------OUTPUT------------
#==============================================================================
print(np.round(Results, decimals=2))

InstalledWind = round(model.Pwind.value,2)
InstalledPV = round(model.Ppv.value, 2)
print("\n----INSTALLED---- \n"
      "Dam Capacity:", Pdam, "MW \n"
      "Wind installiert:", InstalledWind, "MW \n"
      "PV installiert:", InstalledPV, "MW  \n"
      "Für Stromgeneration benötigtes Wasser:", round(sum(model.PowerGeneratingWater[i].value for i in model.T),2), "m^3 \n\n"
      "----COST---- \n"
      "Wind:", "{:0,.2f}".format(model.Cwind * InstalledWind*LifeTimeWind), "€ total |", Cwind, "€ per MW installed \n"
      "PV:", "{:0,.2f}".format(model.Cpv * InstalledPV*LifeTimePv), "€ total |", Cpv, "€ per MW installed \n"
      "Dam", "{:0,.2f}".format(sum((model.PowerGeneratingWater[i].value for i in model.T)) * Pdam * Cdam), "€ total |", Cdam, "€ per MWh used\n"
      "Total:", "{:0,.2f}".format((model.Cwind * InstalledWind * LifeTimeWind + model.Cpv * InstalledPV *LifeTimePv+ sum((model.PowerGeneratingWater[i].value for i in model.T)) * Pdam * Cdam)), "€ \n")

time = datetime.now() - start
print('Problem solved in ' + str(time) + " using " + solver + " as solver.")