
### Display error in calculated electron density
### compared to manufacturers values


import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import math

data = pd.read_csv("relative-e-dens.csv", sep=";")

### PLOT SUBS ON X, RHOS ON Y ###

x = []
for i in range(0,14):
    
    sub = data.iloc[i,0] # Inserts
    x.append(sub)
#------------------
y_calc = []
for i in range(0,14):
    
    rho = float(data.iloc[i,1]) #Calculated Relative Electron Density
    y_calc.append(rho)
#------------------
y_manu = []
for i in range(0,14):
    
    rho = float(data.iloc[i,2]) #Theoretical Relative Electron Density
    y_manu.append(rho)
#------------------

######################## Error-barplot ######################
    
xlen = np.arange(len(x))  # the label locations
width = 0.35  # the width of the bars

rel_err = []
for i in range(0,14):
    
    err = float(data.iloc[i,5]) #Relative error (abs. error / calculated rho)
    rel_err.append(err)

fig, ax = plt.subplots()
rects1 = ax.bar(xlen - width/2, rel_err, width, yerr=0, label='Relative error in %')

# Add some text for labels, title and custom x-axis tick labels, etc.
ax.set_ylabel('Error')
ax.set_title('Calculated \u03C1$_{e}$ error')
ax.set_xticks(xlen)
ax.set_xticklabels(x)
ax.set_axisbelow(True)
ax.yaxis.grid(color='gainsboro')
plt.axis([-1,14,-13,5])
plt.xticks(rotation=30)
plt.yticks(np.arange(-12, 5, 1))
ax.legend()

fig.tight_layout()

plt.show()
