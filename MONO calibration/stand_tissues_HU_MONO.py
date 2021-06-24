### Standard tissues from ICRU Report 44 ###
###      Calculation of reduced HU       ###

import numpy as np
from scipy.optimize import minimize
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

data = pd.read_csv("standard-tissues.csv")

#**# Denominator in A.20 (Goma et al)
wH = 0.112
AH = 1.008
wO = 0.888
AO = 16
k1 = 2.000*10**(-4)
k2 = 3.142*10**(-5)

sum_water = wH/AH*(1+k1+k2) + wO/AO*(8+k1*8**2.86+k2*8**4.62)

calc_HU = []

for i in range(2,64):
    sum_row = 0

    for j in range(1,13):
        wi  = float(data.iloc[i,j])/100          
        Ai  = float(data.iloc[1,j])          
        Zi  = float(data.iloc[0,j])          
        rho = float(data.iloc[i,16]) #relative mass density

        sum_element = sum_element = wi/Ai * (Zi + k1*Zi**2.86 + k2*Zi**4.62)
        sum_row += sum_element

    HU_reduced = rho*(sum_row/sum_water)
    HU = (HU_reduced-1)*1000
    calc_HU.append(HU)

print("*** Calculated HU ***")
for i in range(0,62):
    HUx = calc_HU[i]
    tissue = data.iloc[i+2,0]
    #print (tissue[:13], ":\t", HUx)
    print(HUx)
