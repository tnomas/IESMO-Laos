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
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

#==============================================================================
# --------Initialize-----------
#==============================================================================
#Initialisierung von der 'data.csv', dem Model und der Zeitvariable model.T
#calc_hours = int(input('How many hours should be calculated?: '))
#SpeicherWiederVoll = 1
#SpeicherWiederVoll = int(input('Bis wann soll der Speicher wieder voll sein? (Idealerweise = Calculated Hours): '))
#RiverFlowFactor = float(input('Wieviel Prozent des Flusses von 2017 von der Fluss haben?' ))
Anzeige = int(input("Bis welcher Stunde möchtest du die Werte ansehen: "))
print('Wait...')

data = pd.read_csv('data.csv', header=0, nrows = 8759)
model = ConcreteModel()
model.T = [i for i in range(0,8759)]

#==============================================================================
# ------------Variables------------
#==============================================================================
'''Hier werden die Variablen der einzelnen Technologien und Demands deklariert
To-Do (Wegstreichen wenn fertig):
    Plausible Storage Size von Damm (StorageSize)
    Plausibler Faktor von Wasser zu MWh, aktuell 10,000 m^3 ergeben 1MWh (FactorDam)
    Plausibler Wasserdemand (Einpflegen und als Liste setzen, DemandTotal_Water)
    WaterInflowTotal als Liste und täglicher Einfluss programmieren (WaterInflow)
'''     
    
#Wind
Cwind = 77350 #€/MW/a 
FactorWind = data['windfactor']  #Wind Factor @ hour X

#PV
Cpv = 37120 #€/MWp/a
FactorPv = data['pvfactor']  #Radiation Factor @ hour X

#Dam
StorageSize = 6240000 #m^3 Beschränkt
Pdam = 260 #MW
FactorDam = 20547 #m^3 Water per MW
WaterInflow = data['riverflow_absolut']
Cdam = 0 #€/MWh
DamBalance_Start = StorageSize

#Demand
DemandFactor_Energy = data['energy_demand_normed']
DemandTotal_Energy = 71737.392 #MWh Total

#Costcalculation
DemandEnergy = [(DemandTotal_Energy * DemandFactor_Energy[i]) for i in model.T]
DemandWater = [(data.loc[i,'water_demand_absolut']) for i in model.T]

#==============================================================================
# ------------MODEL CONSTRUCTION------------
#==============================================================================
'''Hier findet die Deklarierung der Variablen und Parameter statt.
Unser Tool optimiert:
    Installierte Leistung Wind
    Installierte Leistung PV
    Stündliche Verwendung des gestauten Wasser zur Energieerzeugung
Simon hat "Sets" nicht gerne gesehen, daher wurde bisher auf diese verzichtet. 
Auch können diese nicht Listen mit doppelten Werten einpflegen: [0,1,2,3,3] ergibt einen Fehler
Es sollten noch alle statischen Variablen die im Modell verwendet werden als Parameter eingepflegt werden (FactorDam usw.)'''

#Parameter
model.Cwind = Param(initialize=Cwind) #Price per MW Wind €/MWh
model.Cpv = Param(initialize=Cpv) #Price per MW PV €/MWh
model.Cdam = Param(initialize=Cdam) #Price of Water €/MWh
model.StorageSize= Param(initialize=StorageSize) #m^3
model.DamBalance_Start = Param(initialize=DamBalance_Start)
model.DamMax = Param(initialize=StorageSize)

#Variablen
model.Pwind = Var(domain=NonNegativeReals) #installed MW Wind
model.Ppv = Var(domain=NonNegativeReals) #installed MW Wind
model.PowerGeneratingWater = Var(model.T, domain=NonNegativeReals) #used M^3 Water for Electricity Generation @ hour X
model.DamBalance = Var(model.T )

#==============================================================================
# ----------Constraint & Objective & Solver------------
#==============================================================================
'''Hier passiert die Magie, Pyomo übernimmt die Constraints und Objective Function'''
#------Objective Function-------
#€/MW * MW + €/MW + m^3(gesamt) * MWh/m^3 * €/MWh, Ziel: Minimieren
def obj_rule(model):
        return(model.Cwind * model.Pwind + model.Cpv * model.Ppv + sum(model.PowerGeneratingWater[i] for i in model.T)/FactorDam * model.Cdam)
    
model.cost = Objective(sense=minimize, rule=obj_rule)
    
#----CONSTRAINTS-------
#Fullfil Demand
#Installierte MW Wind * Faktor Wind + installierte PV * Faktor PV + Benutztes Wasser für Energieerzeugung @ hour X (m^3) * Faktor (MW/m^3) ^>= Verbrauch @ hour X
def DemandEnergy_rule(model, i):
    return (model.Pwind * FactorWind[i] + model.Ppv * FactorPv[i] + model.PowerGeneratingWater[i]/FactorDam >= DemandEnergy[i])

model.EnergyDemand = Constraint(model.T, rule=DemandEnergy_rule)

#Limitation through Turbine Capacity
#Benutztes Wasser für Energieerzeugung @ hour X (m^3) * Faktor (MW/m^3) (= Wasserturbinenleistung @ hour X) <= Maximale Turbinenleistung
def MaxWaterPower_rule(model, i):
    return(model.PowerGeneratingWater[i]/FactorDam <= Pdam)

model.MaxWaterPower = Constraint(model.T, rule=MaxWaterPower_rule)

#Always fullfil Water Demand
'''Hier muss für jede S8tunde immer die aktuelle Situation des Speichers berechnet werden, ein Auslagern auf eine Variable ist wegen Pyomo nicht möglich
darf aber gerne versucht werden! Evtl. muss die Stunde angepasst werden (i+-1)
Wasser für Energieerzeugung <=  Speichergröße (komplette Füllung am Anfang)
                               +Summe bis Stunde X Zufluss
                               -Summe bis Stunde X Wasser für Energieerzeugung
                               -Summe bis Stunde X Wasserverbrauch für Demand'''

def WaterUsage_rule(model, i):
    #return(model.PowerGeneratingWater[i] <= model.StorageSize + sum(WaterInflow[i] for i in range(1,i)) - sum(model.PowerGeneratingWater[i] for i in range(1,i)) - sum(DemandWater[i] for i in range(1,i)))
    if i == 0:
        return(model.DamBalance[i] == model.DamBalance_Start)
    else:
        return(model.DamBalance[i] == model.DamBalance[i-1] + WaterInflow[i] - model.PowerGeneratingWater[i] - DemandWater[i])
        
model.WaterDemand = Constraint(model.T, rule=WaterUsage_rule)

#Never Store more than possible
#Selbes Prinzip wie bei WaterUsage_rule, dieses mal darf es nicht mehr Wasser als die Speichergröße model.StorageSize werden

#model.MaxStorageCapacity = Constraint(model.T, rule=MaxStorageCapacity_rule)

def StorageSameAsBefore_rule (model,i):
    return(model.StorageSize == model.StorageSize + sum(WaterInflow[i] for i in range(1,SpeicherWiederVoll)) - sum(model.PowerGeneratingWater[i] for i in range(1,SpeicherWiederVoll)) - sum(DemandWater[i] for i in range(1,SpeicherWiederVoll)))

#model.StorageSameAsBefore = Constraint(model.T, rule=StorageSameAsBefore_rule)


#------SOLVER---------
opt = SolverFactory('glpk')
model.write('optimization_problem.lp',
         io_options={'symbolic_solver_labels': True})
results = opt.solve(model, tee=True)

#==============================================================================
# ------------Results------------
#==============================================================================
'''Hier wird alles vorbereitet um in dem Dataframe Results zu speichern
Wichtig ist, dass, da die Storagevariable nicht ausgelagert werden kann, Änderungen am Speicher auf in der Variable "Storage" vorgenommen werden'''
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
    
for i in model.T:
    if model.DamBalance[i].value >= model.DamMax.value:
        Storage.append(model.DamMax.value)
    else:
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
Results.index = Results.index + 1 
Results = (np.round(Results, decimals=2))  

#==============================================================================
# ------------Plotting------------
#==============================================================================

ResultsGraph = Results[1:Anzeige]
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
      "Wind:", "{:0,.2f}".format(model.Cwind * InstalledWind), "€ total |", Cwind, "€ per MW installed \n"
      "PV:", "{:0,.2f}".format(model.Cpv * InstalledPV), "€ total |", Cpv, "€ per MW installed \n"
      "Dam", "{:0,.2f}".format(sum((model.PowerGeneratingWater[i].value for i in model.T)) * Pdam * Cdam), "€ total |", Cdam, "€ per MWh used\n"
      "Total:", "{:0,.2f}".format((model.Cwind * InstalledWind + model.Cpv * InstalledPV + sum((model.PowerGeneratingWater[i].value for i in model.T)) * Pdam * Cdam)), "€ \n")

