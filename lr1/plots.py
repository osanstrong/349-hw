import pandas as pd
import matplotlib.pyplot as plt
import thermo as th
import numpy as np

dat = pd.read_csv("boildata.csv")
# print(dat)
PI = 3.141592653589

T = dat["Temperature (C)"]
K = lambda C: C+273.15
dT = T.diff()
dt = dat["Time (sec)"].diff()
print(dT)
print(dt)

# rho = 8933 # kg / m3
COPPER_CASRN = "7440-50-8"
COPPER_M = 63.546 / 1000 # kg / mol
v_Cu_gen = th.VolumeSolid(COPPER_CASRN) # m3 / mol
get_v = lambda T_K: v_Cu_gen(T_K) / COPPER_M # m3 / kg
get_rho = lambda T_K: 1/get_v(T_K) # kg / m3
RHO_RT = get_rho(K(27))
v_RT = get_v(K(27))
R_RT = 4.43 / 100
V_RT = (4/3) * PI * (R_RT**3) # m3

def get_V(T_C):
    return V_RT * get_v(K(T_C)) / v_RT

def get_R(T_C):
    return np.cbrt(get_V(T_C) / (PI * 4/3))
V = V_RT * np.array([get_v(K(T[i]))/v_RT for i in range(len(T))])
R = np.cbrt(V / (PI * 4/3))

k_Cu_gen = th.ThermalConductivitySolid(COPPER_CASRN) # W / m K
k = np.array([k_Cu_gen(K(T[i])) for i in range(len(T))])
MASS = RHO_RT*V_RT
# R = 4.43 / 100 # m
# V_0 = (4/3) * PI * (R**3) # m3
# MASS_0 = 
A = 4 * PI * (R**2) # m2
cp = 385 # J / kg K
cp_Cu_gen = th.HeatCapacitySolid(COPPER_CASRN) # J / molK
get_cp = lambda T_K: cp_Cu_gen(T_K) / COPPER_M # J / kg K
cp = np.array([get_cp(K(T[i])) for i in range(len(T))])
print(f"{cp} {th.HeatCapacitySolid(COPPER_CASRN).units}")
T_SAT = 100

v_Cu_gen = th.VolumeSolid(COPPER_CASRN) # m3 / mol
get_v = lambda T_K: v_Cu_gen(T_K) / COPPER_M # m3 / kg
get_rho = lambda T_K: 1/get_v(T_K) # kg / m3
# get_R = lambda T_K: 

qii = -MASS*cp*(dT/dt)/A
h = qii / (T-T_SAT)
Bi = h * R / k
# print(qii) # W / m2


# Correlations
g = 9.81 # m / s2
sigma = 58.9 / (10**3) # N/m
h_fg = 2257 * 1000 # J / kg
rho_v = 1 / (1.679) # kg / m3
rho_l = 1 / (1.044 / (10**3)) # kg / m3
mu_l = 279 / (10**6) # N s / m2
cp_l = 4217 # J / kg K
cp_v = 2029
Pr_l = 1.76

v_v = (12.02 / (10**6)) / rho_v # m2 / s
k_v = 24.8 / 1000 # W / m K

def hline(qii, label, x = [min(T-T_SAT), max(T-T_SAT)]):
    plt.plot(x, [qii, qii], label=label)
# Free convection

# Nucleate Boiling, eq. 10.5
C_sf = 0.0128
n = 1.0 # Table 10.1
qii_nuc = mu_l * h_fg * np.sqrt(g*(rho_l-rho_v)/sigma) * ((cp_l * (T-T_SAT) / (C_sf*h_fg*(Pr_l**n)))**3)

# Max boiling, eq. 10.6
C = PI / 24 # Zuber constant 
qii_max = C * h_fg * rho_v * (sigma * g  * (rho_l - rho_v) / (rho_v**2))**0.25
# T_max = T[qii.argmax()]
i_max = qii.argmax()
i_crit = abs(qii_nuc - qii_max).argmin()


# Min boiling, eq. 10.7
qii_min = 0.09 * rho_v*h_fg* ((g*sigma*(rho_l-rho_v) / ((rho_l+rho_v)**2))**0.25)

# Film Boiling eq. 10.8, yes for cylinder but we sub C=0.67 for spheres
# D = 2*R
def get_qii_film(T_C):    
    hi_fg = h_fg + 0.8*cp_v*(T_C-T_SAT)
    D = get_R(T_C)

    NuBar_D = 0.67 * (g*(rho_l-rho_v)*hi_fg*(D**3) / (v_v * k_v * (T_C-T_SAT)))**0.25
    hbar_conv = NuBar_D * k_v / get_R(T_C)

    emiss = 0.5 # oxidized copper 
    SB_const = 5.670e-8 # W / m2 K4
    hbar_rad = emiss*SB_const* (K(T_C)**4 - K(T_SAT)**4) / (T_C - T_SAT)

    hbar = hbar_conv + 0.75*hbar_rad
    # hbar = hbar_rad
    qii_film = hbar*(T_C-T_SAT)
    return qii_film

qii_film = np.array([
    get_qii_film(T_C) for T_C in T
])

i_crit2 = abs(qii_film - qii_min).argmin()


print(f"Te of 120: {qii_min:.3f} vs {get_qii_film(100+120)}")


plt.plot(T-T_SAT, qii, label="Observed")
plt.plot(T[i_crit:]-T_SAT, qii_nuc[i_crit:], label="Theoretical nucleate")
# hline(qii_max, "Maximum Convection", x=[10,40])
plt.scatter(T[i_crit]-T_SAT, qii_max, label=f"q''max = {qii_max:.3e} W / m²", c="C1")
# hline(qii_min, "Minimum Convection", x=[100,140])
plt.scatter(T[i_crit2]-T_SAT, qii_min, label=f"q''min = {qii_min:.3e} W / m²", c="C2")
plt.plot(T[:i_crit2]-T_SAT, qii_film[:i_crit2], label="Theoretical film")
plt.xlabel("ΔTe (K)")
plt.ylabel("q'' (W/m²)")
plt.loglog()
plt.legend(loc="best")
plt.show()


# plt.plot(T-T_SAT, Bi)
# plt.xlabel("ΔTe (K)")
# plt.ylabel("Bi")
# plt.loglog()
# plt.legend()
# plt.show()