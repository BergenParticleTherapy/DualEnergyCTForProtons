### PLOT HU vs RSP (ICRU tissues) ###

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import pwlf

data = pd.read_csv("HU_RSP_ICRUt_2_MONO.csv", sep=";")

#########################
### IMPORT DATAPOINTS ###
#########################
HU  = []
RSP = []
##print(data.iloc[0,1])     #For double checking koordinates
##print(data.iloc[1,1])

for i in range(0,58):
    HUx  = float(data.iloc[i,0])
    RSPx = float(data.iloc[i,1])
    HU.append(HUx)
    RSP.append(RSPx)

HUa  = np.array([HU])
RSPa = np.array([RSP])

###########################
### MAKE PIECEWICE PLOT ###
###########################

#Line segment locations
x0 = np.array([min(HU),20,400,max(HU)])

my_pwlf = pwlf.PiecewiseLinFit(HU,RSP)
my_pwlf.fit_with_breaks(x0)

xHat = np.linspace(min(HU), max(HU), num=10000)
yHat = my_pwlf.predict(xHat)

plt.figure()
plt.plot(HU, RSP, 'x')
plt.plot(xHat, yHat, '-')
plt.title('HU vs RSP for ICRU tissues')
plt.xlabel('Calculated HU')
plt.ylabel('Calculated RSP')
plt.show()



#### Print the mathematical equation of the fit ###
from sympy import Symbol
from sympy.utilities import lambdify
x = Symbol('x')


def get_symbolic_eqn(pwlf_, segment_number):
    if pwlf_.degree < 1:
        raise ValueError('Degree must be at least 1')
    if segment_number < 1 or segment_number > pwlf_.n_segments:
        raise ValueError('segment_number not possible')
    # assemble degree = 1 first
    for line in range(segment_number):
        if line == 0:
            my_eqn = pwlf_.beta[0] + (pwlf_.beta[1])*(x-pwlf_.fit_breaks[0])
        else:
            my_eqn += (pwlf_.beta[line+1])*(x-pwlf_.fit_breaks[line])
    # assemble all other degrees
    if pwlf_.degree > 1:
        for k in range(2, pwlf_.degree + 1):
            for line in range(segment_number):
                beta_index = pwlf_.n_segments*(k-1) + line + 1
                my_eqn += (pwlf_.beta[beta_index])*(x-pwlf_.fit_breaks[line])**k
    return my_eqn.simplify()


eqn_list = []
f_list = []
for i in range(my_pwlf.n_segments):
    eqn_list.append(get_symbolic_eqn(my_pwlf, i + 1))
    print('Equation number: ', i + 1)
    print(eqn_list[-1])
    f_list.append(lambdify(x, eqn_list[-1]))

