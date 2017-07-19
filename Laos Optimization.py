#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul 13 17:52:14 2017

@author: nilshoffmann
"""

from pyomo.environ import ConcreteModel, Param, Var, Set, Constraint, Objective, RangeSet
from pyomo.environ import minimize, NonNegativeReals
from pyomo.opt import SolverFactory
import pandas as pd
import numpy as np

#==============================================================================
# --------Initialize-----------
#==============================================================================
data = pd.read_csv('data.csv', index_col = 0, nrows=5)
damdata = pd.read_csv('dam.csv')

model = ConcreteModel()
model.T = data.index
model.DamVariableRange = RangeSet(1,6,1) #damdata.index

#==============================================================================
# ------------Variables------------
#==============================================================================
#Wind
CostWind = 1000000 #€/MW
LifetimeWind = 20 #a
FactorWind = data['windfactor'] #Wind Factor @ hour X

#PV
CostPv = 1000000 #€/MWp
LifetimePv= 20 #a
FactorPv = data['pv_factor']#Radiation Factor @ hour X

#Demand
DemandFactor = data['demand']
DemandTotal = 71737.392

#==============================================================================
# ------------MODEL CONSTRUCTION------------
#==============================================================================
    
#Parameter
model.DemandFactor = Set(initialize=DemandFactor, ordered=True) #Energy demand @ hour X #Possible MWh by Wind @ hour X
model.Cwind = Param(initialize=CostWind/LifetimeWind) #Price per MW Wind
model.Cpv = Param(initialize=CostPv/LifetimePv) #Price per MW PV

#Variablen
model.Pwind = Var(domain=NonNegativeReals) #installed MW Wind
model.Ppv = Var(domain=NonNegativeReals) #installed MW Wind
model.Dam = Var(model.DamVariableRange) #Dam Variable to choose from

#==============================================================================
# ----------Constraint & Objective & Solver------------
#==============================================================================
#------Objective Function-------
        
#Price per MW Wind * installed Wind Capacity + Price per MW Pv * installed PV capacity + Price of Dam installation
def obj_rule(model):
        return(model.Cwind * model.Pwind + model.Cpv * model.Ppv + damdata.loc[model.Dam,'costs'])
    
model.cost = Objective(sense=minimize, rule=obj_rule)


#----CONSTRAINTS-------
#Power of Wind * WindFactor + PV * PVFactor + Waterused @ hour X * Pwater must be bigger than Energy Demand
def DemandEnergy_rule(model, i):
    return (model.Pwind * FactorWind[i] + model.Ppv * FactorPv[i] + damdata.loc[model.Dam,'power'] >= model.DemandFactor[i] * DemandTotal)

model.EnergyDemand = Constraint(model.T, rule=DemandEnergy_rule)


#------SOLVER---------
opt = SolverFactory('glpk')
model.write('optimization_problem.lp',
         io_options={'symbolic_solver_labels': True})
results = opt.solve(model, tee=True)

#==============================================================================
# ------------OUTPUT------------
#==============================================================================
print("")
print("Wind installiert:", round(model.Pwind.value, 2), "MWh")
print("PV installiert:", round(model.Ppv.value, 2), "MWh")
print("Dam installiert:", damdata.loc[model.Dam,'power'])
print("")
for i in model.T:
    print("Stunde", i)
    print("Energy Demand:",model.DemandFactor[i] * DemandTotal, "MWh")
    print("Zusammensetzung:", 
          round(FactorPv[i]*model.Ppv.value,2), "MWh PV |", 
          round(FactorWind[i]*model.Pwind.value,2), "MWh Wind |",
          damdata.loc[model.Dam,'power'], "MWh Dam",)
print("")
print("Wind:", model.Cwind * model.Pwind.value, "€ |", model.Cwind.value, "€ pro MWh" )
print("PV:", model.Cpv * model.Ppv.value, "€ |", model.Cpv.value, "€ pro MWh" )
print("Damm:", damdata.loc[model.Dam,'costs'], "€ |", damdata.loc[model.Dam,'costs']/damdata.loc[model.Dam,'power'], "€ pro MWh" )
print("----")
print("Gesamtkosten:", model.Cwind * model.Pwind.value + model.Cpv * model.Ppv.value)

'''function2(Pdam):
    storage
    
for i in model.T:
     StorageVariable[i] = StorageVariable[i-1] - WasserzuflussFunction(Pdam) + Wasserabflussfunction(Pdam)
    
Wasserzuflussfunction(Pdam):
    if Pdam == 5:
        return(Durchflussgesamt * DatenDämme.Spalte5,Pdam) #km^3
    
Wasserabflussfunction(Pdam):
'''


#Storage
#for i in model.T:
#   
#    
#def WaterUsage_rule (model, i):
#    return(StorageVariable[i] >= model.DemandWater[i])
#
#model.WaterDemand = Constraint(model.T, rule=WaterUsage_rule)