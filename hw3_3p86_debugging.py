import matplotlib.pyplot as plt
import numpy as np

q = 1e5 # W/m³
k_rw = 20 # Wm⁻¹K⁻¹
k_ss = 15
TINF = 25 # C
h = 1000 # W/m²K
r_i = 0.5 # m
r_o = 0.6 # m


r_rw = np.linspace(0, r_i, 1000)

# Fully direct formula for T_rw as requested in 3.86c
def get_T_rw(r):
    return (-q/(6*k_rw)) * r**2   +   q*r_i**2/3 * (0.5/k_rw + (1-r_i/r_o)/k_ss + r_i/(r_o**2*h)) + TINF

T_rw = get_T_rw(r_rw)

C1 = -r_i**3 * q / (3*k_ss)
C2 = TINF + (r_i**3*q/(3*r_o)) * (1/(r_o*h) - 1/k_ss)
C3 = 0
C4 = (q*r_i**2 / 3) * (0.5/k_rw + (1-r_i/r_o)/k_ss + r_i/(r_o**2*h)) + TINF

# T_rw but using the defined coefficients instead of writing out the full formula directly
def get_T_rwb(r):
    return (-q/(6*k_rw)) * r**2 - C3 / r + C4
T_rwb = get_T_rwb(r_rw)


print(f"Temperature at r=0: {get_T_rw(0)}⁰C")
print(f"Temperature at r_i: {get_T_rw(r_i)}⁰C")
print(f"Temperature at r_i: {get_T_rwb(r_i)}⁰C")

r_ss = np.linspace(r_i, r_o, 1000)
T_ss = -C1 / r_ss + C2

plt.plot(r_rw, T_rw, label="T_rw")
# plt.plot(r_rw, T_rwb, label="T_rwb")
plt.plot(r_ss, T_ss, label="T_ss")

plt.xlabel("r (m)")
plt.ylabel("T (⁰C)")
plt.title("Problem 3.86c")
plt.legend()
plt.show()
