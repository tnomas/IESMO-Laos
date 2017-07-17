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
CostPv = 800000 #€/MWp
LifetimePv= 30 #a
FactorPv = [0.2, 0.3, 0.8, 0.4, 0] #Radiation Factor @ hour X

#Water
Storage = 5000 #m^3 in Storage, hourly correction
TurbineLimit = 100 #m^3/h
InflowRain = 40 #m^3 @ hour X
InflowRiver = 200 #m^3 @ hour X

#Demand
DemandEnergy = [5,9000,20,21,300] #MWh @ hour X
DemandWater = 500 #m^3/h @ hour X

''' ------------COST------------ '''
#Costcalculation
Cwind = CostWind / LifetimeWind #Cost/MW Wind, 1 year only
Cpv = CostPv / LifetimePv #Cost/MW PV, 1 year only
Cwater = 10 #€/MWh

''' ------------MODEL CONSTRUCTION------------ '''
model = ConcreteModel()
model.T = Set(initialize = T)

#Parameter
model.DemandWater = Param(initialize=DemandWater) #Water demand @ hour x
model.DemandEnergy = Set(initialize=DemandEnergy, ordered=True) #Energy demand @ hour X
model.FactorWind = Set(initialize=FactorWind, ordered = True)
model.FactorPv = Set(initialize=FactorPv, ordered = True)
model.Cwind = Param(initialize=Cwind) #Possible MWh by Wind @ hour X
model.Cpv = Param(initialize=Cpv) #Possible MWh by PV @ hour X
model.Cwater = Param(initialize=Cwater)

#Variablen
model.Pwind = Var(domain=NonNegativeReals) #installed MW Wind
model.Ppv = Var(domain=NonNegativeReals) #installed MW Wind
model.UsageWater = Var(model.T, bounds = (0,Storage)) #used M^3 Water @ hour X

#Objective Function
def obj_rule(model):
        return(model.Cwind * model.Pwind + model.Cpv * model.Ppv)

model.cost = Objective(sense=minimize, rule=obj_rule)

#Constraints
#MW(installed Wind) * Possible Usage(Wind) + MW(installed PV) + Possible Usage(PV) + 
def demand_rule(model, i):
    return ((model.Pwind * model.FactorWind[i] + model.Ppv * model.FactorPv[i]) >= model.DemandEnergy[i])
    
model.demand = Constraint(model.T, rule=demand_rule)

#Fulfill Drink Water Demand hourly
#Limitation by Water Turbine

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
print("Wind installiert:", model.Pwind.value)
print("PV installiert:", model.Ppv.value)
for i in model.T:
    print("Stunde", i)
    print("Demand:",model.DemandEnergy[i], "MWh")
    print("Zusammensetzung:", model.FactorPv[i]*model.Ppv.value, "aus PV,", model.FactorWind[i]*model.Pwind.value, "aus Wind")