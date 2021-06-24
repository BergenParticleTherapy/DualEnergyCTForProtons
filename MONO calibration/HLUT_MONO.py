### USE HLUT ON MEASURED HU VALUES ###

import numpy as np
import pandas as pd
data = pd.read_csv("Measured_HU_MONO.csv", sep=";")

def HLUT(HU): # HU -> RSP, functions from "plot_HU_RSP_ICRU_2.py"
    if HU < 20:
        RSP = 0.00109025235442672*HU + 1.1055812541744
    elif 20 < HU < 400:
        RSP = 0.000827605725033786*HU + 1.0793165912351
    elif HU > 400:
        RSP = 0.000594775916103616*HU + 1.12588255302114
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

##HU = float(data.iloc[1,0])
##print(HU)
