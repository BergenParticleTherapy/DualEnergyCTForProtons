#####              USES  SECT               #####
### Calculates HU from gammex using k1 and k2 ###
#####                                       #####

import numpy as np
from scipy.optimize import minimize
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

data = pd.read_csv("datablad_gammex.csv", sep = ";")

#**# Denominator in A.20 (Goma et al)
wH = 0.112
AH = 1.008
wO = 0.888
AO = 16               ######################################
k1 = 2.00000*10**(-4) ### K-VALUES FROM SECT CALIBRATION ###
k2 = 2.04346*10**(-5) ######################################

sum_water = wH/AH*(1+k1+k2) + wO/AO*(8+k1*8**2.86+k2*8**4.62)
    #**#
    
sum_tot = 0

#print(data.iloc[3,5])
for i in range (2,18):          # itterate through rows (CB2-30% - Breast)
    sum_row = 0
        
    for j in range (1,10):      # itterate through columns (H - Ca)
        wi  = float(data.iloc[i,j])          
        Ai  = float(data.iloc[1,j])          
        Zi  = float(data.iloc[0,j])          
        rho = float(data.iloc[i,11]) #relative mass density         

        sum_element = wi/Ai * (Zi + k1*Zi**2.86 + k2*Zi**4.62)
        sum_row += sum_element
            
    calc_HU_red = rho*(sum_row/sum_water)
    HU = (calc_HU_red -1)*1000
    print(HU)
