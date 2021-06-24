import numpy as np
import matplotlib.pyplot as plt
from sklearn.linear_model import LinearRegression


rcs_l = np.array([1.239, 1.181, 1.535, 1.184, 1.393]).reshape((-1,1))
rcs_h = np.array([1.082, 1.106, 1.177, 1.062, 1.128])


model = LinearRegression().fit(rcs_l,rcs_h)

r_sq = model.score(rcs_l,rcs_h)
print("Coefficient of determination: ", r_sq)

print("Intercept: ", model.intercept_)
print("Slope: ", model.coef_)
