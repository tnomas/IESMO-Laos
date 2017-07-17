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

#Time
T = [1,2,3,4,5]

''' ------------Variables------------ '''
#Wind
CostWind = 1000000 #€/MW
LifetimeWind = 20 #a
FactorWind = [0.5, 0.3, 0.1, 0.6, 0.2] #Wind Factor @ hour X

#PV
CostPv = 400000 #€/MWp
LifetimePv= 30 #a
FactorPv = [0.2, 0.3, 0.8, 0.4, 0.02] #Radiation Factor @ hour X

#Water
Storage = 10000000 #m^3 in Storage, hourly correction
Pwater = 0.001 #MW/m^3
TurbineLimit = 480 #MW
InflowRain = 40 #m^3 @ hour X
InflowRiver = 200 #m^3 @ hour X

#Demand
DemandEnergy = [500, 400, 300, 450, 310] #MWh @ hour X
DemandWater = [200, 300, 400, 301, 100] #m^3/h @ hour X

''' ------------COST------------ '''
#Costcalculation
Cwind = CostWind #Cost/MW Wind, 1 year only
Cpv = CostPv #Cost/MW PV, 1 year only
Cwater = 1 #€/MWh

''' ------------MODEL CONSTRUCTION------------ '''
model = ConcreteModel()
model.T = Set(initialize = T)

#Parameter
#model.DemandWater = Set(initialize=DemandWater, ordered=True) #Water demand @ hour x
model.DemandEnergy = Set(initialize=DemandEnergy, ordered=True) #Energy demand @ hour X
model.FactorWind = Set(initialize=FactorWind, ordered = True)
model.FactorPv = Set(initialize=FactorPv, ordered = True)
model.Cwind = Param(initialize=Cwind) #Possible MWh by Wind @ hour X
model.Cpv = Param(initialize=Cpv) #Possible MWh by PV @ hour X
#model.Cwater = Param(initialize=Cwater)
#model.Storage = Param(initialize=Storage, mutable = True)
model.TurbineLimit = Param(initialize=TurbineLimit)
model.Pwater = Param(initialize=Pwater)

#Variablen
model.Pwind = Var(domain=NonNegativeReals) #installed MW Wind
model.Ppv = Var(domain=NonNegativeReals) #installed MW Wind
model.UsageWater = Var(model.T, domain=NonNegativeReals) #used M^3 Water @ hour X

''' ----------Constraint & Objective ------------ '''
#------Objective Function-------
def obj_rule(model):
        return(model.Cwind * model.Pwind + model.Cpv * model.Ppv + sum(model.UsageWater[i] for i in model.T) * Pwater * Cwater)
model.cost = Objective(sense=minimize, rule=obj_rule)


#----CONSTRAINTS-------
#Power of Wind * WindFactor + PV * PVFactor + Waterused @ hour X * Pwater must be bigger than Energy Demand
def DemandEnergy_rule(model, i):
    return (model.Ppv * model.FactorPv[i] + model.Pwind * model.FactorWind[i] + model.UsageWater[i] * Pwater >= model.DemandEnergy[i])
model.Energy = Constraint(model.T, rule=DemandEnergy_rule)


#Limitation through Water Turbine
def MaxWaterPower_rule(model, i):
    return(model.UsageWater[i] * model.Pwater  <= model.TurbineLimit)
model.MaxWaterPower = Constraint(model.T, rule=MaxWaterPower_rule)

#Solver
opt = SolverFactory('glpk')
model.write('optimization_problem.lp',
         io_options={'symbolic_solver_labels': True})

''' ------------OUTPUT------------ '''
#Ausgabe
#model.pprint()

results = opt.solve(model, tee=True)

#print(results)
print("")
print("Wind installiert:", round(model.Pwind.value, 2), "MWh")
print("PV installiert:", round(model.Ppv.value, 2), "MWh")
print("Gesamt benötigtes Wasser:", round(sum(model.UsageWater[i].value for i in model.T),2), "m^3")
print("")
for i in model.T:
    print("Stunde", i)
    print("Demand:",model.DemandEnergy[i], "MWh")
    print("Zusammensetzung:", 
          round(model.FactorPv[i]*model.Ppv.value,2), "MWh PV |", 
          round(model.FactorWind[i]*model.Pwind.value,2), "MWh Wind |", 
          round(model.UsageWater[i].value*Pwater,2), "MWh Wasser")
    print(round(model.UsageWater[i].value,2), "m^3 gebrauchtes Wasser")
print("")
print("Wind:", model.Cwind * model.Pwind.value, "€ |", Cwind, "€ pro MWh" )
print("PV:", model.Cpv * model.Ppv.value, "€ |", Cpv, "€ pro MWh" )
print("Wasser:", sum((model.UsageWater[i].value for i in model.T)) * Pwater * Cwater, "€ |", Cwater, "€ pro MWh" )
print("----")
print("Gesamtkosten:", model.Cwind * model.Pwind.value + model.Cpv * model.Ppv.value + sum((model.UsageWater[i].value for i in model.T)) * Pwater * Cwater, "€")