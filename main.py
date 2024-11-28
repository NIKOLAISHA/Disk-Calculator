
import calfromdatabase as dt
# from scipy.integrate import simpson
from numpy import pi
from numpy import zeros
from matplotlib.pyplot import plot, xlabel, ylabel, title, legend, grid

# from numpy import array

# read file

fp = open("experdata.txt","r",encoding="utf-8")
lines = fp.readlines()
numlines = len(lines)
fp.close()


mu0 = 0.5 #Coefficient Poisson without PLASTIC DEFORMATION

def trIntegral(Y,x):
    numPoints = len(Y)
    result = zeros(numPoints)
    for i in range(1, numPoints):
        result[i] = result[i-1] + (Y[i] + Y[i-1]) / 2 * (x[i] - x[i-1])
    return result

for i in range(0,numlines-1):
    if lines[i] == "Boundary conditions\n":
        #           Boundary conditions
        templist=lines[i+3].split()
        n = float(templist[0])
        SigmaRa = float(templist[1])*98100      #Convert from kgf/cm^2 to Pa
        SigmaRb = float(templist[2])*98100      #Convert from kgf/cm^2 to Pa

        # rouw2 = float(templist[3])
        tau = float(templist[4])  
    if lines[i] == "Points\n":
        numPoints = int(lines[i+1].split()[0])
        r = zeros(numPoints)
        h = zeros(numPoints)
        T = zeros(numPoints)
        for j in range(i+4,i+4+numPoints):
            templist = lines[j].split()
            r[j-i-4] = float(templist[0])/100
            h[j-i-4] = float(templist[1])/100
            T[j-i-4] = float(templist[2])
            
rouw2 = 8480 * (n * pi / 30) ** 2  #Density of material EI698 IS 8480 kg/m3
# rouw2 = 1015
#calculate parameters module Young and others

E0 = zeros(numPoints)
Alpha = zeros(numPoints)
SigmaDl = zeros(numPoints)
SigmaEl = zeros(numPoints)
SigmaEk = zeros(numPoints)
SigmaR = zeros(numPoints)
SigmaTau = zeros(numPoints)
psi = zeros(numPoints)

for i in range(0,numPoints):
    E0[i] = dt.getE(T[i])*100000000000   #Convert to Pa
    Alpha[i] = dt.getAlpha(T[i])*0.000001
    SigmaDl[i] = dt.getSigmaDl(T[i])*1000000000
    SigmaEl[i] = dt.getSigmaEl(T[i])*1000000000
    print(str(T[i])+"\t"+str(E0[i])+"\t"+str(Alpha[i])+"\t"+str(SigmaDl[i])+"\t"+str(SigmaEl[i]))

#calculate strength in irretions

err = 100 #Initialization
psi_ = zeros(numPoints)
psi = zeros(numPoints)
epsilonP = zeros(numPoints)
# tempComp = rea

for irr in range(0,10000):
    
    if irr >= 1:
        tempComp = (SigmaEl / SigmaEk - 1 + (1 - SigmaEl/SigmaEk)\
                    /(SigmaDl - SigmaEl) * (E0 * 0.15 - SigmaEl)) * SigmaEk / E0
    # a
    for i in range(0,numPoints-1):
        if SigmaEk[i] <= SigmaEl[i]:
            pass
        
        elif SigmaEk[i] > SigmaEl[i] and SigmaEk[i] < SigmaDl[i]:
            if tempComp[i] > epsilonP[i]:
                epsilonP[i] = tempComp[i] 
             
        elif SigmaEk[i] > SigmaDl[i]:
            print("Out of range")
    
    if irr >= 1:
        psi_ = epsilonP / ( SigmaEk / E0 ) 
    
    mu = mu0            #mu is constant

    psi = psi_          #psi equal to psi_
    E = E0 / (1 + psi)
    G = E / 2 / (1 + mu)
    
    #Get K
    K = zeros(numPoints)
    integral0 = trIntegral(r ** mu * Alpha * T, r)
    K = (1 + mu) * E / r ** (1 + mu) * integral0 - E * Alpha * T + \
        r[0] ** (1 + mu) * Alpha[0] * T[0] * E / (r ** (1 + mu))
    A = zeros(numPoints)    
    
    integral1 = trIntegral( E * h / r ** (2 + mu), r)
#    print(simpson(E * h / r ** (2 + mu),x=r))
    integral2 = trIntegral(h * K / r  , r)
    integral3 = trIntegral(r * h, r)
    integral4 = trIntegral( r ** mu * SigmaR / E, r )
    integral5 = trIntegral(h / r * (((1 - mu ** 2) * E / r ** (1 + mu) * integral4 - (1 - mu) * SigmaR)) , r )  
    # integral6 = trIntegral(h / r * ((1 - mu ** 2) * E / r ** (1 + mu) * (integral4[numPoints-1] - (1 - mu) * SigmaR)) , r )
    
    A = 1 / integral1[numPoints-1] * (SigmaRb * h[numPoints - 1] - SigmaRa \
        * h[0] - integral2[numPoints-1] + rouw2 * integral3[numPoints-1] \
            - integral5[numPoints-1])
    
    # A = r[0] ** (1 + mu) / E[0] * (SigmaTau[0] - mu * SigmaR[0])
    
    tempSigmaR = 1 / h * (integral2 - rouw2 * integral3 + integral5 + h[0] * \
            SigmaRa + A * trIntegral(E * h / r ** (2 + mu) , r))    
    tempSigmaTau = mu * tempSigmaR + K + (1 - mu ** 2) * E / r ** (1 + mu) * trIntegral(r ** mu * tempSigmaR / E , r)\
        + A * E / r ** (1 + mu)
    tempSigmaEk = (tempSigmaR ** 2 + tempSigmaTau ** 2 - tempSigmaR * tempSigmaTau) ** 0.5
    
    err = max(abs(tempSigmaEk - SigmaEk))
    
    SigmaR = tempSigmaR
    SigmaTau = tempSigmaTau
    SigmaEk = tempSigmaEk
    
    print(irr , err)
    if err < 0.0000001:
        break
        
plot(r * 1000 , SigmaR / 1000000, color = "r", marker = "*")    
plot(r * 1000 , SigmaTau / 1000000, color = "b",  marker = "+") 
plot(r * 1000 , SigmaEk / 1000000, color = "k",  marker = "^") 
xlabel("Radius, mm")
ylabel("Strength, MPa")
legend(['Radial Strength','Tangential Strength','Eqvivalent Strength'])
grid()
title("Strenth Plot")
