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

# ==============================================================================
# --------Initialize-----------
# ==============================================================================

WaterFlowFactor = 10
data = pd.read_csv('data.csv', header=0)
model = ConcreteModel()
model.T = [i for i in range(0, 8760)]

# ==============================================================================
# ------------Variables------------
# ==============================================================================

# Wind
Cwind = 77350  # €/MW/a
FactorWind = data['windfactor']  # Wind Factor @ hour X
LifeTimeWind = 20

# PV
Cpv = 37120  # €/MWp/a
FactorPv = data['pvfactor']  # Radiation Factor @ hour X
LifeTimePv = 25

# Dam
StorageSize = 6240000  # m^3
Pdam = 260  # MW
FactorDam = 20547  # m^3 Water per MW
WaterInflow = 0.01 * WaterFlowFactor * data['riverflow_absolut']
Cdam = 0.0001  # €/MWh
DamBalance_Start = StorageSize / 2

# Demand
DemandFactor_Energy = data['energy_demand_normed']
DemandTotal_Energy = 71737.392  # MWh Total

# Costcalculation
DemandEnergy = [(DemandTotal_Energy * DemandFactor_Energy[i]) for i in model.T]
DemandWater = data['water_demand_absolut']

print('Wait...')

# ==============================================================================
# ------------MODEL CONSTRUCTION------------
# ==============================================================================
# Parameter
model.Cwind = Param(initialize=Cwind)  # Price per MW Wind €/MWh
model.Cpv = Param(initialize=Cpv)  # Price per MW PV €/MWh
model.Cdam = Param(initialize=Cdam)  # Price of Water €/MWh
model.StorageSize = Param(initialize=StorageSize)  # m^3
model.DamBalance_Start = Param(initialize=DamBalance_Start)

# Variablen
model.Pwind = Var(domain=NonNegativeReals)  # installed MW Wind
model.Ppv = Var(domain=NonNegativeReals)  # installed MW Wind
model.PowerGeneratingWater = Var(model.T, domain=NonNegativeReals)  # used M^3 Water for Electricity Generation @ hour X
model.DamBalance = Var(model.T, bounds=(20000, StorageSize))

# ==============================================================================
# ----------Constraint & Objective & Solver------------
# ==============================================================================

# ------Objective Function-------
# €/MW * MW + €/MW + m^3(gesamt) * MWh/m^3 * €/MWh, Ziel: Minimieren
def obj_rule(model):
        return(model.Cwind * model.Pwind + model.Cpv * model.Ppv + 
               sum(model.PowerGeneratingWater[i] for i in model.T)/FactorDam *
               model.Cdam)

model.cost = Objective(sense=minimize, rule=obj_rule)

# ----CONSTRAINTS-------
# Fullfil Demand
def DemandEnergy_rule(model, i):
    return (model.Pwind * FactorWind[i] + model.Ppv * FactorPv[i] +
            model.PowerGeneratingWater[i]/FactorDam >= DemandEnergy[i])

model.EnergyDemand = Constraint(model.T, rule=DemandEnergy_rule)


# Limitation through Turbine Capacity
def MaxWaterPower_rule(model, i):
    return(model.PowerGeneratingWater[i]/FactorDam <= Pdam)

model.MaxWaterPower = Constraint(model.T, rule=MaxWaterPower_rule)


# Calculate Storage
def WaterUsage_rule(model, i):
    if i == 0:
        return(model.DamBalance[i] == model.DamBalance_Start)
    else:
        return(model.DamBalance[i] == model.DamBalance[i-1] + WaterInflow[i] -
               model.PowerGeneratingWater[i] - DemandWater[i])

model.WaterDemand = Constraint(model.T, rule=WaterUsage_rule)

# Dont use more Water than in Storage
def GenerationWaterLimit_rule(model, i):
    if i == 0:
        return(model.PowerGeneratingWater[i] <= model.DamBalance_Start)
    else:
        return(model.PowerGeneratingWater[i] <= model.DamBalance[i])

model.GenerationWaterLimit = Constraint(model.T, rule=GenerationWaterLimit_rule)


# Same Storage in the End level like @ hour 1
def DamEndLvl_rule(model, i):
    return(model.DamBalance[8759] == model.DamBalance_Start)

model.DamEndLvl = Constraint(model.T, rule=DamEndLvl_rule)

# ------SOLVER---------
opt = SolverFactory('glpk')
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
for i in model.T:
    Hour.append(i)
    WindOutput.append(FactorWind[i]*model.Pwind.value)
    PVOutput.append(FactorPv[i]*model.Ppv.value)
    DamOutput.append(model.PowerGeneratingWater[i].value/FactorDam)
    GenerationWaterUsed.append(model.PowerGeneratingWater[i].value)
    WaterBalance.append(WaterInflow[i] - model.PowerGeneratingWater[i].value - DemandWater[i])
    River.append(WaterInflow[i])
    Storage.append(model.DamBalance[i].value)

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
                        "Water Balance": pd.Series(WaterBalance)},
                        columns=['Hour','Energy Demand','Wind Output','PV Output','Dam Output','Factor Wind','Factor PV','Storage','Turbine Water','Water Demand','River','Water Balance'])
                        
Results = (np.round(Results, decimals=2))  

#==============================================================================
# ------------Plotting------------
#==============================================================================

AnzeigeAnfang = 0 #int(input("Ab welcher Stunde möchtest du die Werte ansehen: "))
AnzeigeEnde = 8760 #int(input("Bis welcher Stunde möchtest du die Werte ansehen: "))

ResultsGraph = Results[AnzeigeAnfang:AnzeigeEnde]
demand = ResultsGraph['Energy Demand']
y1 = ResultsGraph['PV Output']
y2 = y1 + ResultsGraph['Dam Output']
y3 = y2 + ResultsGraph['Wind Output']
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

#==============================================================================
# ------------OUTPUT------------
#==============================================================================
print(np.round(Results, decimals=2))

InstalledWind = round(model.Pwind.value,2)
InstalledPV = round(model.Ppv.value, 2)
print("\n----INSTALLED---- \n"
      "Dam Capacity:", Pdam, "MW \n"
      "Wind installiert:", InstalledWind, "MWh \n"
      "PV installiert:", InstalledPV, "MWh  \n"
      "Für Stromgeneration benötigtes Wasser:", round(sum(model.PowerGeneratingWater[i].value for i in model.T),2), "m^3 \n\n"
      "----COST---- \n"
      "Wind:", "{:0,.2f}".format(model.Cwind * InstalledWind*LifeTimeWind), "€ total |", Cwind, "€ per MW installed \n"
      "PV:", "{:0,.2f}".format(model.Cpv * InstalledPV*LifeTimePv), "€ total |", Cpv, "€ per MW installed \n"
      "Dam", "{:0,.2f}".format(sum((model.PowerGeneratingWater[i].value for i in model.T)) * Pdam * Cdam), "€ total |", Cdam, "€ per MWh used\n"
      "Total:", "{:0,.2f}".format((model.Cwind * InstalledWind * LifeTimeWind + model.Cpv * InstalledPV *LifeTimePv+ sum((model.PowerGeneratingWater[i].value for i in model.T)) * Pdam * Cdam)), "€ \n")
