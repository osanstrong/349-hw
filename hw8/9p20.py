import thermo as th
from thermo import HeatCapacityLiquid, ViscosityLiquid, VolumeLiquid, ThermalConductivityLiquid
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

print(f"Water viscosity: {visc(T):.10f}")
print(f"Water viscosity: {7*10**-7:.10f} (manual)")

def Ra(L, Ti, Ts):
    Ta = (Ti+Ts)/2
    return (g*beta(Ta)*abs(Ti-Ts)*L**3)/(visc(T)*alpha(T))

L = 0.2 # m
print(f"Ra for {K(50)} water and {K(30)} surface: {Ra(L, K(50), K(30))}")
