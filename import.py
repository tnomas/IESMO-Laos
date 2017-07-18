#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul 18 16:03:02 2017

@author: thomasbreitenstein
"""
from pyomo.environ \
    import ConcreteModel, Param, Var, Constraint, Objective, Set, RangeSet, NonNegativeReals
from pyomo.environ import maximize, minimize
from pyomo.opt import SolverFactory
import pandas as pd

#==============================================================================
# ACHTUNG: nur hierfür erstellt. Einmal löschen beim kopieren
#==============================================================================
model = ConcreteModel(name='Districtheat')


#==============================================================================
# jetzt wirds interessant ;-) 
#==============================================================================

data = pd.read_csv('data.csv', index_col = 0, nrows=24)
damdata = pd.read_csv('dam.csv')

model.T = data.index
model.demand = Param(initialize=data['demand'])

#==============================================================================
# und so weiter und so fort... 
#==============================================================================
