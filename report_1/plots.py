import matplotlib.pyplot as plt
import numpy as np
import json

input = json.load(open("input.json", 'r'))
# Copper Material Properties
rho = 8933 # kg / m^3
c_p = 385 # J / kg*K
if "k_Wm-1K-1" in input:
    k = input["k_Wm-1K-1"] # W / m*K, default of 401
else:
    k = 401
alph = k / (rho*c_p) # (m^2 / s)

# Fluid properties
h = input["h_Wm-2K-1"] # W/m^2*K
T_INF = input["T_INF_C"] # C

# Experiment Setup
R = input["R_m"] # m
T_0 = input["T_0_C"] # C
t_final = input["total_time_s"] # s
V = 4/3 * np.pi * R**3
A = 4 * np.pi * R**2


# analytic plot A1
t = np.linspace(0, t_final, 1000)
tau = rho*c_p*V / (h*A)

T = T_INF + (T_0-T_INF)*np.exp(-t / tau)


# print(f"T({t[-1]}): {T[-1]}")

# plt.plot(t, T, label=f"Analytic T(t), D={2*R}m, h={h:.3f}W/(m^2K)")
# plt.xlabel("t (s)")
# plt.ylabel("T (degrees C)")
# plt.legend()
# plt.show()

#### Plot A2 ####
# Base
tau = rho*c_p*V / (h*A)
T_A1 = T_INF + (T_0-T_INF)*np.exp(-t / tau)
# Double A
A2 = A*2
tau = rho*c_p*V / (h*A2)
T_2A = T_INF + (T_0-T_INF)*np.exp(-t / tau)
# Double R (halve A/V ratio)
R2 = 2*R
V8 = 4/3 * np.pi * R2**3
A4 = 4 * np.pi * R2**2
tau = rho*c_p*V8 / (h*A4)
T_2R = T_INF + (T_0-T_INF)*np.exp(-t / tau)
# plt.plot(t, T_A1, label=f"Base dimensions: D={2*R}m, A={A*1e6:.2f}cm^2")
# plt.plot(t, T_2A, label=f"Base radius, double A: D={2*R}m, A={A2*1e6:.2f}cm^2")
# plt.plot(t, T_2R, label=f"Double radius: D={2*R2}m, A={A4*1e6:.2f}cm^2")
# plt.xlabel("t (s)")
# plt.ylabel("T (degrees C)")
# plt.legend()
# plt.show()


# Double h
h2 = h*2
tau = rho*c_p*V / (h2*A)
T_h2 = T_INF + (T_0-T_INF)*np.exp(-t / tau)

plt.plot(t, T_A1, label=f"Base h ({h:.3f}W/(m^2K))")
plt.plot(t, T_2A, label=f"Double h ({h2:.3f}W/(m^2K))")
plt.xlabel("t (s)")
plt.ylabel("T (degrees C)")
plt.legend()
plt.show()


# Plot analytic vs different dt of FD

# Plot analytic vs different dr

