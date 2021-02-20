#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
This is a simple example for LP formulation with Gurobi-Python
This problem is called Resource Assignment Problem (RAP)

This code is entirely based on the following link:
    https://www.gurobi.com/resource/tutorial-mixed-integer-linear-programming/
"""
# import gurobi library
from gurobipy import *

# resources and job sets [list of resources and jobs]
R = ['Carlos', 'Joe', 'Monika']
J = ['Tester', 'JavaDeveloper', 'Architect']

# matching score data and cost in a multidict
"""
Multidict allows you to initialize one or more dictionaries in 
a single statement. Here, ('R','J') is the key for dictionary and
ms indicates the matching score for each combination
C indicates the cost of each combination
"""
combinations, ms, C = multidict({
    ('Carlos', 'Tester'): [53, 1],
    ('Carlos', 'JavaDeveloper'): [27, 1],
    ('Carlos', 'Architect'): [13, 1],
    ('Joe', 'Tester'): [80, 2],
    ('Joe', 'JavaDeveloper'): [47, 2],
    ('Joe', 'Architect'): [67, 2],
    ('Monika', 'Tester'): [53, 3],
    ('Monika', 'JavaDeveloper'): [73, 3],
    ('Monika', 'Architect'): [47, 3]
    })

B = 5 # budget limit

# Inirialize model
m = Model('RAP') # RAP is the name of this model

# Create decision variables for the RAP model
x = m.addVars(combinations, vtype=GRB.BINARY, name='assign')

# Create gap variables to show that we have the option "not assign a job"
g = m.addVars(J, name='gap')


# Create jobs constraints
jobs = m.addConstrs((x.sum('*', j) + g[j] == 1 for j in J), 'job')
"""
Considering: 
    Carlos=1; Joe=2; Monika=3
    Tester=1; JavaDeveloper=2; JavaDeveloper=3
The above line creates the following constraints:
    x(1,1) + x(1,2) + x(1,3) + g1 == 1
    x(2,1) + x(2,2) + x(2,3) + g2 == 1
    x(3,1) + x(3,2) + x(3,3) + g3 == 1
"""
# Create resources constraints
resources = m.addConstrs((x.sum(r, '*') <= 1 for r in R), 'resource')

# Create budget constrain
budget = m.addConstr((x.prod(C) <= B), 'budget')

# Define big M for penalizing the option of "not assign a job"
bigM = 100

# Define the objective function to maximize the total matching score
m.setObjective(x.prod(ms) - bigM*g.sum(), GRB.MAXIMIZE)

# Save model for inspection
m.write('RAP.lp')

# Run optimization engine
m.optimize()

# Display optimal values of decision variables
for i in m.getVars():
    if (abs(i.x) > 1e-5):
        print(i.varName, i.x)
        
# Display the optimized value of the objective function
print('optimal objective function', m.objVal)

# Display total matching scores and total costs after optimization
total_matching_scores = 0
costs = 0
for [r, j] in combinations:
    if (abs(x[r, j].x) > 1e-5):
        total_matching_scores = total_matching_scores + ms[r, j]*x[r, j].x
        costs = costs + C[r, j]*x[r, j].x
        
print('total matching scores', total_matching_scores)
print('total costs', costs)