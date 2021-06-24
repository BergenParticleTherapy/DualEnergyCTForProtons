### USE HLUT ON MEASURED HU VALUES ###

import numpy as np
import pandas as pd
data = pd.read_csv("Measured_HU.csv", sep=";")

def HLUT(HU): # HU -> RSP, functions from "plot_HU_RSP_ICRU_2.py"
    if HU < 20:
        RSP = 0.00107776081155286*HU + 1.09800661310044
    elif 20 < HU < 400:
        RSP = 0.00087734136329856*HU + 1.07796466827501
    elif HU > 400:
        RSP = 0.000733057372973769*HU + 1.10682146633997
    return RSP

###Get HU values to be converted to RSP


HU_list = []
for i in range(0,16):
    HU = float(data.iloc[i,0])
    HU_list.append(HU)
print(HU_list)


#Convert to RSP
RSP_list = []
for i in range(0,16):
    HU = float(data.iloc[i,0])
    RSP = HLUT(HU)
    RSP_list.append(RSP)
    print(RSP)


