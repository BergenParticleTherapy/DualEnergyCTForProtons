
### Display DIFF between measured HU and calculated HU
### MONO


import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import math

data = pd.read_csv("HUt_vs_HUc.csv")
print(data)
### PLOT SUBS ON X, HU-diff ON Y ###

#Get label names for x-axis
x = []
for i in range(0,17):
    
    sub = data.iloc[i,0] # Inserts
    x.append(sub)

###################### Error-barplot ######################
    
xlen = np.arange(len(x))  # the label locations
width = 0.35  # the width of the bars

rel_err = []
for i in range(0,17):
    
    err = float(data.iloc[i,8]) #Relative error (abs. error / calculated rho)
    rel_err.append(err)

fig, ax = plt.subplots()
rects1 = ax.bar(xlen - width/2, rel_err, width, yerr=0, label='Diff')

# Add some text for labels, title and custom x-axis tick labels, etc.
ax.set_ylabel('measured HU - calculated HU')
ax.set_title('Calculated HU error (MONO)')
ax.set_xticks(xlen)
ax.set_xticklabels(x)
ax.set_axisbelow(True)
ax.yaxis.grid(color='gainsboro')
plt.axis([-1,15,-13,5])
plt.xticks(rotation=30)
plt.yticks(np.arange(-12, 85, 4))
#ax.legend()

fig.tight_layout()

plt.show()
