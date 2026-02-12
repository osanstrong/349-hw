import matplotlib.pyplot as plt
import numpy as np


# Constants
# T1 = 261 #C
# T2 = 211 #C
T1 = 830 #C
T2 = 361 #C
TINF = 25 #C
LA = 0.03
LB = 0.03 #M
LC = 0.02

kA = 25 # W/mK
# hA = 1000 # W/m²K
hA = 0
kC = 50 
hC = 1000 # W/m²K



q = 4e6 #W/m³
k = 15.35 #W/mK

# C1 = -833.33 #K/m
C1 = (T2-T1)/(2*LB)
C2 = 0.5*T1 + 0.5*T2 + (0.5*q/k)*(LB**2) #K

# C3 = (T1-TINF) / (kA/hA + LA)
C3 = 0
C4 = T1 + C3*LB

C5 = (TINF - T2) / (LC + kC/hC)
C6 = T2 - C5*LB

print(f"C2 equal to {C2}K")


xB = np.linspace(-LB, LB, 1000)
TBx = (-0.5*q/k) * xB*xB   +   C1*xB  + C2

xA = np.linspace(-LB-LA, -LB, 1000)
TAx = C3*xA + C4

xC = np.linspace(LB, LB+LC, 1000)
TCx = C5*xC + C6

plt.plot(xA, TAx, label="TA")
plt.plot(xB, TBx, label="TB")
plt.plot(xC, TCx, label="TC")

plt.xlabel("x (m)")
plt.ylabel("T (C)")
plt.legend()

# plt.title("Problem 3.68b")
plt.title("Problem 3.68c")
plt.show()