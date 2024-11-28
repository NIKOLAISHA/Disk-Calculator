# -*- coding: utf-8 -*-
"""
Created on Mon Nov 25 09:48:37 2024

@author: 18059
"""

#calculate from material.txt

fp = open("material.txt","r",encoding="utf-8")

lines = fp.readlines()
numLines = len(lines)

ET = [] 
E = []
Alpha = []
AlphaT = []
SigmaDl = []
SigmaDlT = []
SigmaEl = []
SigmaElT = []

#READ MATERIAL FILE

for i in range(0,numLines-1):
    if lines[i] == "E\n":
        j=int(lines[i+1].split()[0])
        for k in range(i+2,i+2+j):
            ET.append(float(lines[k].split()[0])-273.15) #convert string to float
            E.append(float(lines[k].split()[1]))
    elif lines[i] == "ALPHA\n":
        j=int(lines[i+1].split()[0])
        for k in range(i+2,i+2+j):
            AlphaT.append(float(lines[k].split()[0])-273.15)
            Alpha.append(float(lines[k].split()[1]))
    elif lines[i] == "SIGMADL\n":
        j=int(lines[i+1].split()[0])
        for k in range(i+2,i+2+j):
            SigmaDlT.append(float(lines[k].split()[0])-273.15)
            SigmaDl.append(float(lines[k].split()[1]))
    elif lines[i] == "SIGMAEL\n":
        j=int(lines[i+1].split()[0])
        for k in range(i+2,i+2+j):
            SigmaElT.append(float(lines[k].split()[0])-273.15)
            SigmaEl.append(float(lines[k].split()[1]))
        
fp.close()
    


def getE(T):
    for i in range(0,len(ET)-1):
        if T > ET[i] and T<= ET[i+1]:
            E0 = (E[i+1]-E[i])/(ET[i+1]-ET[i])*(T-ET[i])+E[i]
            return E0
    return "Out Of Range"

def getAlpha(T):
    for i in range(0,len(ET)-1):
        if T > AlphaT[i] and T<= AlphaT[i+1]:
            E0 = (Alpha[i+1]-Alpha[i])/(AlphaT[i+1]-AlphaT[i])*(T-AlphaT[i])+Alpha[i]
            return E0
    return "Out Of Range"

def getSigmaDl(T):
    for i in range(0,len(ET)-1):
        if T > SigmaDlT[i] and T<= SigmaDlT[i+1]:
            E0 = (SigmaDl[i+1]-SigmaDl[i])/(SigmaDlT[i+1]-SigmaDlT[i])*(T-SigmaDlT[i])+SigmaDl[i]
            return E0
    return "Out Of Range"

def getSigmaEl(T):
    for i in range(0,len(ET)-1):
        if T > SigmaElT[i] and T<= SigmaElT[i+1]:
            E0 = (SigmaEl[i+1]-SigmaEl[i])/(SigmaElT[i+1]-SigmaElT[i])*(T-SigmaElT[i])+SigmaEl[i]
            return E0
    return "Out Of Range"