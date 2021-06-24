# Gradient decent minimization to determine k1 and k2

import numpy as np
from scipy.optimize import minimize
from scipy.stats import linregress
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

data = pd.read_csv("datablad_gammex_MONO.csv", sep=";")
measured = pd.read_csv("Measured_HU_MONO.csv", sep=";")

def objective(k):
    k1 = k[0]
    k2 = k[1]

    x_ = list()
    y_ = list()

    #**# Denominator in A.20 (Goma et al)
    wH = float(data.iloc[16,1])
    AH = float(data.iloc[1,1])
    wO = float(data.iloc[16,4])
    AO = float(data.iloc[1,4])

    sum_water = wH/AH*(1+k1+k2) + wO/AO*(8+k1*8**2.86+k2*8**4.62)
    #**#
    
    sum_tot = 0

    for i in range (2,18):          # itterate through rows (CB2-30% - Breast)
        sum_row = 0
        
        for j in range (1,10):      # itterate through columns (H - Ca)
            wi  = float(data.iloc[i,j])          
            Ai  = float(data.iloc[1,j])          
            Zi  = float(data.iloc[0,j])          
            rho = float(data.iloc[i,11]) #relative mass density         

            sum_element = wi/Ai * (Zi + k1*Zi**2.86 + k2*Zi**4.62)
            sum_row += sum_element
            
        calc_HU = rho*(sum_row/sum_water)
        meas_HU = float(measured.iloc[i-2,1])

        if kDRAW:
            x_.append(meas_HU)
            y_.append(calc_HU)

        diff_squared = (calc_HU - meas_HU)**2

        sum_tot += diff_squared

    if kDRAW:
        fig = plt.figure(figsize=(10,5))
        plt.subplot(121)
        plt.scatter(x_, y_, label="data")
        xlin = np.linspace(min(x_), max(x_), 100)
        plt.xlabel("Measured HU")
        plt.ylabel("Calculated HU")
        plt.ylim([0, 2.5])
        plt.title(f"k1 = {k1:.3e}; k2 = {k2:.3e}; error = {sum_tot:.3e}")
        plt.plot(xlin, xlin, "r", label=f"x=y")
        plt.legend()
        
        plt.subplot(122)
        plt.scatter(x_, [a-b for a,b in zip(x_,y_)], label="calculated - measured")
        xlin = np.linspace(min(x_), max(x_), 100)
        plt.xlabel("Measured HU")
        plt.ylabel("Calculated - measured HU")
        plt.ylim([-0.1, 0.1])
        plt.show()
    
    return sum_tot

# Initial gues-values
k0 = [5.3*10**(-4),2.3*10**(-5)]
kDRAW = True
print(objective(k0))

kDRAW = False
# Bounds (Goma et al, p. 5)
b1 = (2*10**(-4),6*10**(-3))
b2 = (3*10**(-6),6*10**(-4))
bnds = (b1,b2)

# Solve:
sol = minimize(objective, k0, method='SLSQP', bounds=bnds)
print(sol)

kDRAW = True
print(objective(sol.x))


