
### Calulate relative electron densities rho_e
### Calculate I-values for tissue subs


import pandas as pd
import numpy as np
from scipy.constants import c as c
import matplotlib.pyplot as plt
import math

data = pd.read_csv("datablad_gammex.csv", sep=';')
#print("*** Datablad_Gammex 467 ***")
#print(data)

### Calculation of relative electron densities
NA = 6.0221409*10**23   # Avrogados

wH = float(data.iloc[16,1])
AH = float(data.iloc[1,1])
ZH = float(data.iloc[0,1])
wO = float(data.iloc[16,4])
AO = float(data.iloc[1,4])
ZO = float(data.iloc[0,4])

n_ew = 1*NA*(wH*(ZH/AH) + wO*(ZO/AO))   # Electron density of water
rho_list = []

print("***Electron densities***")
for i in range (2,18):                  # rows (CB2-30% - Breast)
    rho_m = float(data.iloc[i,11])      # mass densities
    sum = 0
    
    for j in range (1,10):              # columns (H - Ca)
        wi  = float(data.iloc[i,j])          
        Ai  = float(data.iloc[1,j])          
        Zi  = float(data.iloc[0,j])
        sum += wi*(Zi/Ai)               # the sum in A.3 Goma
    
    n_e = rho_m*NA*sum
    rho_e = n_e/n_ew                    # Relative electron densities

    substitute = data.iloc[i,0]
    print (substitute[:5],":\t",rho_e)
    rho_list.append(rho_e)              ### Save rho_e's to list

#print(rho_list)

### Calculation of mean I-values for all subs
#print("***Mean I-values***")
e = math.e
I_list = []
for i in range (2,18):                  # rows (CB2-30% - Breast)
    sum1 = 0
    sum2 = 0
    
    for j in range (1,10):              # columns (H - Ca)
        lni = math.log(float(data.iloc[18,j]))
        wi  = float(data.iloc[i,j])          
        Ai  = float(data.iloc[1,j])          
        Zi  = float(data.iloc[0,j])
        sum1 += wi*(Zi/Ai)*lni
        sum2 += wi*(Zi/Ai)

    lni_mean = (sum1/sum2)
    I_list.append(lni_mean)             ### Save I-values to list
    
##print("***I-values***")
##print(I_list)
#print(len(I_list))

####### Calculation of RSPs
MeV = 1e6
Iw = 78# eV
Ek = 100*MeV             # 100 MeV
me = 0.511*MeV           # electron mass  ###OBS! Burde være i eV?
mp = 938*MeV                         ###OBS! Burde være i eV?
b  = 0.4282              # Velocity of the 100 MeV proton [c]
#print (b)

S_ew = math.log(2*me*c**2*b**2)-math.log(1-b**2)-b**2-math.log(Iw) #Denomin

print("***RSP***")
for i in range (0,16):
    Ix = I_list[i]
    pe = rho_list[i]
    substitute = data.iloc[i+2,0]
    S_e = math.log(2*me*c**2*b**2)-math.log(1-b**2)-b**2-math.log(Ix)

    RSP = pe*(S_e/S_ew)
    #print (substitute[:5],":\t",RSP) # og lagre disse et sted. now we gottem
    #print(RSP)
