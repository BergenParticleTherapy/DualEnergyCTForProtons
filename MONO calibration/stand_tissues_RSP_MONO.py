### Standard tissues from ICRU Report 44 ###
###          Calculation of RSP          ###

import pandas as pd
import numpy as np
from scipy.constants import c as c
import matplotlib.pyplot as plt
import math

data = pd.read_csv("standard-tissues.csv")
#print(data)

### Reading of relative electron densities
rho_list = []

for i in range(2,64):                   #Rows (Adipose# - Water)
    rho_tissue = float(data.iloc[i,16]) #Columns (16 = rho_e_rel)
    rho_list.append(rho_tissue)
   


### Calculation of mean I-values for all tissues
e = math.e
I_list = []
for i in range (2,64):
    sum1 = 0
    sum2 = 0

    for j in range (1,13):              #Columns (H - I)
        lni = math.log(float(data.iloc[64,j]))
        wi  = (float(data.iloc[i,j]))/100          
        Ai  = float(data.iloc[1,j])          
        Zi  = float(data.iloc[0,j])
        sum1 += wi*(Zi/Ai)*lni
        sum2 += wi*(Zi/Ai)

    lni_mean = (sum1/sum2)
    I_list.append(lni_mean)
##print("***I-values***")


######### Calculation of RSPs #########
MeV = 1e6
Iw = 78# eV
Ek = 100*MeV             # 100 MeV
me = 0.511*MeV           
mp = 938*MeV            
b  = 0.4282              # Velocity of the 100 MeV proton [c]

S_ew = math.log(2*me*c**2*b**2)-math.log(1-b**2)-b**2-math.log(Iw) #Denomin

print("***RSP***")
RSP_list = []
for i in range (0,62):
    Ix = I_list[i]
    pe = rho_list[i]
    substitute = data.iloc[i+2,0]
    S_e = math.log(2*me*c**2*b**2)-math.log(1-b**2)-b**2-math.log(Ix)

    RSP = pe*(S_e/S_ew)
    RSP_list.append(RSP)
    #print (substitute[:13],":\t",RSP)
    print(RSP)
    


