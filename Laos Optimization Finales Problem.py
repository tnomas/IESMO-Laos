#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''
Created on Thu Jul 13 17:52:14 2017
@author: nilshoffmann
'''

from pyomo.environ import ConcreteModel, Param, Var, Set, Constraint, Objective
from pyomo.environ import minimize, NonNegativeReals
from pyomo.opt import SolverFactory
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from datetime import datetime
import time

#==============================================================================
# --------Initialize-----------
#==============================================================================
#Initialisierung von der 'data.csv', dem Model und der Zeitvariable model.T
calc_hours = int(input('How many hours should be calculated? '))
print('Wait...')
start = datetime.now()


data = pd.read_csv('data.csv', header=0, index_col = 0, nrows = calc_hours)
model = ConcreteModel()
model.T = [i for i in data.index]

#==============================================================================
# ------------Variables------------
#==============================================================================
'''Hier werden die Variablen der einzelnen Technologien und Demands deklariert
To-Do (Wegstreichen wenn fertig):
    Plausible Wind, PV und Dammkosten (Dammkosten von 0?, CostWind, CostPv)
    Plausible Lebenszeit von PV und Wind (LifetimeWind, LifetimePv)
    Plausible Storage Size von Damm (StorageSize)
    Plausibler Faktor von Wasser zu MWh, aktuell 10,000 m^3 ergeben 1MWh (FactorDam)
    Plausible Leistung des Dams (Pdam)
    Plausibler Wasserdemand (Einpflegen und als Liste setzen, DemandTotal_Water)
    WaterInflowTotal als Liste und täglicher Einfluss programmieren (WaterInflow)
'''     
    
#Wind
CostWind = 1000000 #€/MW Muss noch realisitisch gemacht werden
LifetimeWind = 20 #a
FactorWind = [0]+[(data.loc[i,'windfactor']*8760) for i in model.T]  #Wind Factor @ hour X

#PV
CostPv = 10000 #€/MWp
LifetimePv= 20 #a
FactorPv = [0]+[(data.loc[i,'pv_factor']*8760) for i in model.T]  #Radiation Factor @ hour X

#Dam
StorageSize = 90000 #m^3 Beschränkt
StorageVolume = 50000
Pdam = 100 #MW
FactorDam = 10000 #m^3 Water per MW
FactorWaterInflow = [0] + [data.loc[i,'riverflow_normed'] for i in model.T]
WaterInflowTotal = 1529431.498 * 1000 #m^3
Cdam = 0 #€/MWh

#Demand
DemandFactor_Energy = data['demand']
DemandTotal_Energy = 71737.392 #MWh Total
DemandFactor_Water = data['water_demand']
DemandTotal_Water =  45376800    #m^3 Total anually

#Costcalculation
Cwind = CostWind/LifetimeWind #€/MW Wind, 1 year only
Cpv = CostPv/LifetimePv #€/MW PV, 1 year only
DemandEnergy = [0]+[(DemandTotal_Energy * DemandFactor_Energy[i]) for i in data.index]
DemandWater = [0]+[(DemandTotal_Water * DemandFactor_Water[i]) for i in data.index]

#Excessive Demand
CExcess = 9999999

#==============================================================================
# ------------MODEL CONSTRUCTION------------
#==============================================================================
'''Hier findet die Deklarierung der Variablen und Parameter statt.
Unser Tool optimiert:
    Installierte Leistung Wind
    Installierte Leistung PV
    Stündliche Verwendung des gestauten Wasser zur Energieerzeugung

Simon hat 'Sets' nicht gerne gesehen, daher wurde bisher auf diese verzichtet. 
Auch können diese nicht Listen mit doppelten Werten einpflegen: [0,1,2,3,3] ergibt einen Fehler
Es sollten noch alle statischen Variablen die im Modell verwendet werden als Parameter eingepflegt werden (FactorDam usw.)'''

#Parameter
model.Cwind = Param(initialize=Cwind) #Price per MW Wind €/MWh
model.Cpv = Param(initialize=Cpv) #Price per MW PV €/MWh
model.Cdam = Param(initialize=Cdam) #Price of Water €/MWh
model.StorageSize= Param(initialize=StorageSize) #m^3
model.CExcess = Param(initialize=CExcess)


#Variablen
model.Pwind = Var(domain=NonNegativeReals) #installed MW Wind
model.Ppv = Var(domain=NonNegativeReals) #installed MW Wind
model.PowerGeneratingWater = Var(model.T, domain=NonNegativeReals) #used M^3 Water for Electricity Generation @ hour X
model.PExcess = Var(model.T, domain=NonNegativeReals)
model.StorageVolume = Var(model.T, domain=NonNegativeReals, bounds=(0, StorageSize))
model.StorageExcess = Var(model.T)

#==============================================================================
# ----------Constraint & Objective & Solver------------
#==============================================================================
'''Hier passiert die Magie, Pyomo übernimmt die Constraints und Objective Function'''
#------Objective Function-------
#€/MW * MW + €/MW + m^3(gesamt) * MWh/m^3 * €/MWh, Ziel: Minimieren
def obj_rule(model):
        return(model.Cwind * model.Pwind + model.Cpv * model.Ppv + sum(model.PowerGeneratingWater[i] for i in model.T)/FactorDam * Cdam \
                                                                       + model.CExcess * (sum(model.PExcess[i] for i in model.T)+ sum(model.StorageExcess[i] for i in model.T)))
    
model.cost = Objective(sense=minimize, rule=obj_rule)
    
#----CONSTRAINTS-------
#Fullfil Demand
#Installierte MW Wind * Faktor Wind + installierte PV * Faktor PV + Benutztes Wasser für Energieerzeugung @ hour X (m^3) * Faktor (MW/m^3) ^>= Verbrauch @ hour X
def DemandEnergy_rule(model, i):
    return (model.Pwind * FactorWind[i] + model.Ppv * FactorPv[i] + model.PowerGeneratingWater[i]/FactorDam + model.PExcess[i] >= DemandEnergy[i])

model.EnergyDemand = Constraint(model.T, rule=DemandEnergy_rule)

#Limitation through Turbine Capacity
#Benutztes Wasser für Energieerzeugung @ hour X (m^3) * Faktor (MW/m^3) (= Wasserturbinenleistung @ hour X) <= Maximale Turbinenleistung
def MaxWaterPower_rule(model, i):
    return(model.PowerGeneratingWater[i]/FactorDam <= Pdam)

model.MaxWaterPower = Constraint(model.T, rule=MaxWaterPower_rule)

#Always fullfil Water Demand
'''Hier muss für jede Stunde immer die aktuelle Situation des Speichers berechnet werden, ein Auslagern auf eine Variable ist wegen Pyomo nicht möglich
darf aber gerne versucht werden! Evtl. muss die Stunde angepasst werden (i+-1)
Wasser für Energieerzeugung <=  Speichergröße (komplette Füllung am Anfang)
                               +Summe bis Stunde X Zufluss
                               -Summe bis Stunde X Wasser für Energieerzeugung
                               -Summe bis Stunde X Wasserverbrauch für Demand'''

def WaterUsage_rule(model, i):
    return(model.PowerGeneratingWater[i] <= model.StorageSize + sum((FactorWaterInflow[i] * WaterInflowTotal) for i in range(1,i)) - sum(model.PowerGeneratingWater[i] for i in range(1,i)) - sum(DemandWater[i] for i in range(1,i)))

model.WaterDemand = Constraint(model.T, rule=WaterUsage_rule)

#Never Store more than possible
#Selbes Prinzip wie bei WaterUsage_rule, dieses mal darf es nicht mehr Wasser als die Speichergröße model.StorageSize werden
def MaxStorageCapacity_rule (model, i):
    if 1 < i < 8760:
        return(model.StorageVolume[i] == model.StorageVolume[i-1] + FactorWaterInflow[i] * WaterInflowTotal - model.PowerGeneratingWater[i] - DemandWater[i] + model.StorageExcess[i])
        #return(model.StorageSize >= model.StorageSize + sum((FactorWaterInflow[i] * WaterInflowTotal) for i in range(1,i)) - sum(model.PowerGeneratingWater[i] for i in range(1,i)) - sum(DemandWater[i] for i in range(1,i)))
    else:
        return(model.StorageVolume[i] == StorageVolume + FactorWaterInflow[i] * WaterInflowTotal - model.PowerGeneratingWater[i] - DemandWater[i]+ model.StorageExcess[i])

model.MaxStorageCapacity = Constraint(model.T, rule=MaxStorageCapacity_rule)

#------SOLVER---------
opt = SolverFactory('gurobi')
model.write('optimization_problem.lp',
         io_options={'symbolic_solver_labels': True})
results = opt.solve(model, tee=True)

#==============================================================================
# ------------Results------------
#==============================================================================
'''Hier wird alles vorbereitet um in dem Dataframe Results zu speichern
Wichtig ist, dass, da die Storagevariable nicht ausgelagert werden kann, Änderungen am Speicher auf in der Variable 'Storage' vorgenommen werden'''
Hour = []
WindOutput = []
EnergyDemand = []; WaterDemand = []
PVOutput = []; DamOutput = []
GenerationWaterUsed = []; Storage = [0]
WaterInflow = [0]
PExcess = []
StorageVolume = []
for i in model.T:
    Hour.append(i)
    WindOutput.append(round(FactorWind[i]*model.Pwind.value,2))
    PVOutput.append(round(FactorPv[i]*model.Ppv.value,2))
    DamOutput.append(round(model.PowerGeneratingWater[i].value/FactorDam,2))
    GenerationWaterUsed.append(round(model.PowerGeneratingWater[i].value,2))
    Storage.append(round((model.StorageSize + sum((FactorWaterInflow[i] * WaterInflowTotal) for i in range(1,i)) - sum(model.PowerGeneratingWater[i].value for i in range(1,i)) -sum(DemandWater[i] for i in range(1,i))),2))
    WaterInflow.append(round((FactorWaterInflow[i] * WaterInflowTotal),2))
    PExcess.append(round((model.PExcess[i].value),2))
    StorageVolume.append(round((model.StorageVolume[i].value),2))

Results = pd.DataFrame({'Hour': pd.Series(Hour), 
                        'Energy Demand': pd.Series(DemandEnergy[1:]),
                        'Water Demand': pd.Series(DemandWater[1:]),
                        'Wind Output': pd.Series(WindOutput),
                        'PV Output': pd.Series(PVOutput),
                        'Dam Output': pd.Series(DamOutput),
                        'Turbine Water': pd.Series(GenerationWaterUsed),
                        'Factor Wind': pd.Series(FactorWind[1:]),
                        'Factor PV': pd.Series(FactorPv[1:]),
                        'Storage': pd.Series(Storage[1:]),
                        'Water Inflow': pd.Series(WaterInflow[1:]),
                        'Excess': pd.Series(PExcess),
                        'Storage Level': pd.Series(StorageVolume)},
                        columns=['Hour','Energy Demand','Water Demand','Wind Output','PV Output','Dam Output','Turbine Water','Factor Wind','Factor PV','Storage','Water Inflow','Excess','Storage Level'])
Results.index = Results.index + 1


Results['Storage'] = Results['Storage'] / 10000 #nötig? Sonst sind die Werte gigantisch groß, da Wasser und nicht Leistung


#==============================================================================
# ------------Plotting------------
#==============================================================================

#Results.plot.area()
ResultsGraph = Results[1:25]
demand = ResultsGraph['Energy Demand']
y1 = ResultsGraph['PV Output']
y2 = y1 + ResultsGraph['Dam Output']
y3 = y2 + ResultsGraph['Storage']
y4 = y3 + ResultsGraph['Wind Output']
y5 = y4 + ResultsGraph['Excess']

fig = plt.figure()
ax = fig.add_subplot(111)
ax.fill_between(ResultsGraph.Hour, 0, y1, label = 'PV', color = 'm', alpha=.7)
ax.fill_between(ResultsGraph.Hour, y1, y2, label = 'Dam', color = 'g', alpha=.7)
ax.fill_between(ResultsGraph.Hour, y2, y3, label = 'Storage', color = 'orange', alpha=.7)
ax.fill_between(ResultsGraph.Hour, y3, y4, label = 'Wind', color = 'b', alpha=.7)
ax.fill_between(ResultsGraph.Hour, y4, y5, label = 'Excess', color = 'r', alpha=.7)
ax.plot(ResultsGraph.Hour, demand, label = 'Energy Demand', color = 'black')
ax.yaxis.set_label_position('right')
ax.set_ylabel('MW')
ax.yaxis.set_label_position('left')

box = ax.get_position()
ax.set_position([box.x0, box.y0 + box.height * 0.1, box.width, box.height * 0.9])

plt.legend(fontsize = 'small', loc='upper center', bbox_to_anchor=(0.5, -0.15), fancybox=True, ncol=5)
plt.xlabel('Hours')
plt.title('Energy output')
plt.savefig('graphical_output/energy.pdf', dpi=150)


#==============================================================================
# ------------OUTPUT------------
#==============================================================================
InstalledWind = round(model.Pwind.value,2)
InstalledPV = round(model.Ppv.value, 2)
print('\n----INSTALLED---- \n'
      'Dam Capacity:', Pdam, 'MW \n'
      'Wind installiert:', InstalledWind, 'MWh \n'
      'PV installiert:', InstalledPV, 'MWh  \n'
      'Für Stromgeneration benötigtes Wasser:', round(sum(model.PowerGeneratingWater[i].value for i in model.T),2), 'm^3 \n\n'
      '----COST---- \n'
      'Wind:', '{:0,.2f}'.format(model.Cwind * InstalledWind), '€ total |', Cwind, '€ per MW installed \n'
      'PV:', '{:0,.2f}'.format(model.Cpv * InstalledPV), '€ total |', Cpv, '€ per MW installed \n'
      'Dam', '{:0,.2f}'.format(sum((model.PowerGeneratingWater[i].value for i in model.T)) * Pdam * Cdam), '€ total |', Cdam, '€ per MWh used\n'
      'Total:', '{:0,.2f}'.format((model.Cwind * InstalledWind + model.Cpv * InstalledPV + sum((model.PowerGeneratingWater[i].value for i in model.T)) * Pdam * Cdam)), '€ \n')
print(Results)

time = datetime.now() - start
print('Runtime: ' + str(time) + " for a calculation of " + str(calc_hours) + " hours.")

#plt.show()

