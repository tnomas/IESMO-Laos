#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul 13 17:52:14 2017

@author: nilshoffmann
"""

from pyomo.environ import ConcreteModel, Param, Var, Set, Constraint, Objective
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
DamVariable = [1,2,3,4,5] #damdata.index

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
FactorPv = [0.5, 0.3, 0.8, 0.4, 0.0] #Radiation Factor @ hour X

#

#Demand
DemandFactor = data['demand']
DemandTotal = 71737.392
#[400, 500, 300, 450, 310] #MWh @ hour X

#==============================================================================
# ------------MODEL CONSTRUCTION------------
#==============================================================================
#Parameter
model.DemandFactor = Set(initialize=DemandFactor, ordered=True) #Energy demand @ hour X
model.FactorWind = Set(initialize=FactorWind, ordered = True) #Possible MWh by Wind @ hour X
model.FactorPv = Set(initialize=FactorPv, ordered = True) #Possible MWh by PV @ hour X
model.Cwind = Param(initialize=CostWind/LifetimeWind) #Price per MW Wind
model.Cpv = Param(initialize=CostPv/LifetimePv) #Price per MW PV

#Variablen
model.Pwind = Var(domain=NonNegativeReals) #installed MW Wind
model.Ppv = Var(domain=NonNegativeReals) #installed MW Wind
model.Dam = Var(bounds=(1,2)) #Dam Variable to choose from

#==============================================================================
# ----------Constraint & Objective & Solver------------
#==============================================================================
#------Objective Function-------

def pricedam(index):
#    if (0 < index <= 2):
#        return(10)
    #Oder das:
    return(damdata.loc[index,'costs'])
        
#Price per MW Wind * installed Wind Capacity + Price per MW Pv * installed PV capacity + Price of Dam installation
def obj_rule(model):
        return(model.Cwind * model.Pwind + model.Cpv * model.Ppv + pricedam(model.Dam))
    
model.cost = Objective(sense=minimize, rule=obj_rule)


#----CONSTRAINTS-------
#Power of Wind * WindFactor + PV * PVFactor + Waterused @ hour X * Pwater must be bigger than Energy Demand
def DemandEnergy_rule(model, i):
    return (model.Pwind * model.FactorWind[i] + model.Ppv * model.FactorPv[i] >= model.DemandFactor[i] * DemandTotal)

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
print("")
for i in model.T:
    print("Stunde", i)
    print("Energy Demand:",model.DemandFactor[i] * DemandTotal, "MWh")
    print("Zusammensetzung:", 
          round(model.FactorPv[i]*model.Ppv.value,2), "MWh PV |", 
          round(model.FactorWind[i]*model.Pwind.value,2), "MWh Wind |")
print("")
print("Wind:", model.Cwind * model.Pwind.value, "€ |", model.Cwind, "€ pro MWh" )
print("PV:", model.Cpv * model.Ppv.value, "€ |", model.Cpv, "€ pro MWh" )
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