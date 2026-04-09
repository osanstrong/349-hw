import thermo as th
from thermo import HeatCapacityLiquid, ViscosityLiquid, VolumeLiquid, ThermalConductivityLiquid
import matplotlib.pyplot as plt

T = 310 # K?
def K(C):
    return C+273.15
P = 101325 #Pa
M = 18.015e-3 # kg / mol
g = 9.81 # m/s^2
# P = 1.0 # Atm
CASRN_water="7732-18-5"
mu = ViscosityLiquid(CASRN=CASRN_water)
vol = VolumeLiquid(CASRN=CASRN_water)

# print(f"Units: {vol.units}")
cp = HeatCapacityLiquid(CASRN=CASRN_water)
print(f"cp units: {cp.units}")
k = ThermalConductivityLiquid(CASRN=CASRN_water)
def beta(T):
    return (1/vol(T, P)) * vol.T_dependent_property_derivative(T)
print(f"Water β: {beta(T)}")
print(f"Water β: {361.9/1e6:.10f}")
def alpha(T):
    # Thermal diffusivity
    return k(T, P)*vol(T, P)/cp(T)
def rho(T):
    return M / vol(T,P)
print(f"Water rho: {rho(T)} kg/m³")
def visc(T):
    return mu(T, P)/rho(T)
def Pr(T):
    return visc(T)/alpha(T)

print(f"Water viscosity: {visc(T):.10f}")
print(f"Water viscosity: {7*10**-7:.10f} (manual)")

def Ra(L, Ti, Ts):
    Ta = (Ti+Ts)/2
    return (g*beta(Ta)*abs(Ti-Ts)*L**3)/(visc(T)*alpha(T))

# General case Nu, Eq. 9.26
def Nu(Ra, Pr):
    return (0.825 + (0.387*(Ra**(1/6)))/((1 + (0.492/Pr)**(9/16))**(8/27)))**2

def Nu_turb(Ra):
    C = 0.1
    n = 1/3
    return C*(Ra**n)

def calc_h(L, Ti, Ts):
    Ra_val = Ra(L, Ti, Ts)
    Ta = (Ti+Ts)/2
    Nu_val = Nu(Ra_val, Pr(T))
    h = Nu_val * k(Ta, P) / L
    return h

L = 0.2 # m
Ra_val = Ra(L, K(50), K(30))
print(f"Ra for {K(50)} water and {K(30)} surface: {Ra_val:,}")
Nu_val = Nu(Ra_val, Pr(T))
# Nu_val = Nu_turb(Ra_val)
print(f"Pr = {Pr(T)}")
print(f"-> Nu = {Nu_val}")
print(f"  -> h = {Nu_val * k(K(40), P) / L} W/m²K")


TH = K(50)
TC = K(10)

def calc_overall_h(TH, TC, debug=True):
    TS_init = (TH+TC)/2
    TS = TS_init
    hH = calc_h(L, TH, TS)
    hC = calc_h(L, TC, TS)
    if debug: print(f"0  TS: {TS:.5f}, hH={hH:.3f}, hC={hC:.3f}")

    for i in range(10):
        hH = calc_h(L, TH, TS)
        hC = calc_h(L, TC, TS)
        TS = (hH*TH + hC*TC) / (hC + hH)
        if debug: print(f"{i}  TS: {TS:.5f}, hH={hH:.3f}, hC={hC:.3f}")

    h_tot = hC * (TS-TC) / (TH-TC)
    if debug: print(f"Total h: {h_tot} W/m²K")
    if debug: print(f"q = {hH*(TH-TS)}\n  = {hC*(TS-TC)}\n  = {h_tot*(TH-TC)}")
    return h_tot

h_tot_9p20a = calc_overall_h(TH, TC)

TH_vals_C = [i for i in range(20, 61)]
h_tot_vals = [calc_overall_h(K(TH_val_C), TC, debug=False) for TH_val_C in TH_vals_C]
# print(h_tot_vals)

plt.plot(TH_vals_C, h_tot_vals, label="Overall h")
plt.xlabel("Hot process temperature (⁰C)")
plt.ylabel("Overall h (W/m²K)")
plt.legend()
plt.show()
