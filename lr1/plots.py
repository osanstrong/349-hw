import pandas as pd
import matplotlib.pyplot as plt

dat = pd.read_csv("boildata.csv")
# print(dat)

T = dat["Temperature (C)"]
dT = T.diff()
dt = dat["Time (sec)"].diff()
print(dT)
print(dt)

rho = 8933 # kg / m3
R = 4.43 / 100 # m
PI = 3.141592653589
V = (4/3) * PI * (R**3) # m3
A = 4 * PI * (R**2) # m2
cp = 385 # J / kg K
T_SAT = 100

qii = -rho*V*cp*(dT/dt)/A
print(qii) # W / m2
plt.plot(T-T_SAT, qii)
plt.xlabel("ΔTe (K)")
plt.ylabel("q'' (W/m²)")
plt.legend()
plt.show()