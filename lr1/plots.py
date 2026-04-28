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

v_l = (279e-6) / rho_l
v_v = (12.02 / (10**6)) / rho_v # m2 / s
k_l = 680 / 1000
k_v = 24.8 / 1000 # W / m K

def hline(qii, label, x = [min(T-T_SAT), max(T-T_SAT)]):
    plt.plot(x, [qii, qii], label=label)


# Free convection
Beta = 750.1e-6 # 1 / K
def get_Gr(T_C):
    return g*Beta*(T_C-T_SAT)*((2*R_RT)**3) / (v_l**2)
def get_Ra(T_C):
    return get_Gr(T_C)*Pr_l
def get_NuBarFree(T_C):
    return 2 + (0.589*(get_Ra(T_C)**0.25)) / ((1 + (0.469/Pr_l)**(9/16))**(4/9))
def get_hbarFree(T_C):
    Nu = get_NuBarFree(T_C)
    return Nu * k_l / (2*R_RT)
def get_qiiFree(T_C):
    return get_hbarFree(T_C)*(T_C-T_SAT)


T_free = np.linspace(min(T), 105, 30)
qii_free = get_qiiFree(T_free)


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
# T_nuc = np.linspace


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


fig, ax1 = plt.subplots()




ax1.plot(T-T_SAT, qii, label="Observed")
ax1.plot(T[i_crit:]-T_SAT, qii_nuc[i_crit:], label="Theoretical nucleate")
# hline(qii_max, "Maximum Convection", x=[10,40])
ax1.scatter(T[i_crit]-T_SAT, qii_max, label=f"q''max = {qii_max:.3e} W / m²", c="C1")
# hline(qii_min, "Minimum Convection", x=[100,140])
ax1.scatter(T[i_crit2]-T_SAT, qii_min, label=f"q''min = {qii_min:.3e} W / m²", c="C2")
ax1.plot(T[:i_crit2]-T_SAT, qii_film[:i_crit2], label="Theoretical film")
ax1.plot(T_free-T_SAT, qii_free, label="Theoretical Free Convection")
ax1.set_xlabel("ΔTe (K)")
ax1.set_ylabel("q'' (W/m²)")


ax1.loglog()

qmin, qmax = ax1.get_ylim()
Te_valid = 240
ax1.plot([Te_valid, Te_valid], [qmin, qmax], label="Bi=0.1")
ax1.set_ylim(qmin, qmax)


# ax1.set_zorder(1)
ax1.legend()


plt.show()
fig, ax2 = plt.subplots()

# ax2 = plt.twinx(ax1)
ax2.plot(T-T_SAT, Bi, label="Biot Number")
ax2.plot([min(T-T_SAT), max(T-T_SAT)], [0.1, 0.1], label="Te ~ 240 K")

# ax1.loglog()
ax2.loglog()
bmin, bmax = ax2.get_ylim()
ax2.plot([Te_valid, Te_valid], [bmin, bmax], label="Bi=0.1")
ax2.set_ylim(bmin, bmax)
ax2.legend(loc="best")
ax2.set_xlabel("ΔTe (K)")
ax2.set_ylabel("Bi (-)")

plt.show()


# plt.plot(T-T_SAT, Bi)
# plt.xlabel("ΔTe (K)")
# plt.ylabel("Bi")
# plt.loglog()
# plt.legend()
# plt.show()