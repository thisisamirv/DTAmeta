#!/usr/local/bin/python3

import sympy as sym
import pandas as pd

data = pd.read_csv('Data.csv', sep=",")

TPs = []
TNs = []
FPs = []
FNs = []

for i in enumerate(data['SEN']):
    print(i[0])
    TP,TN,FP,FN = sym.symbols('TP,TN,FP,FN')
    eq1 = sym.Eq(((TP+TN)/(TP+TN+FP+FN)), (data['ACC'][i[0]])/100)
    eq2 = sym.Eq((TP/(TP+FN)), (data['SEN'][i[0]])/100)
    eq3 = sym.Eq((TN/(TN+FP)), (data['SPE'][i[0]])/100)
    eq4 = sym.Eq((TP+TN+FP+FN), data['N'][i[0]])
    result = sym.solve([eq1,eq2, eq3, eq4],(TP,TN,FP,FN))
    TPs.append(result[TP])
    TNs.append(result[TN])
    FPs.append(result[FP])
    FNs.append(result[FN])
