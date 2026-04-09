import numpy as np
from thermo import HeatCapacityGas, ViscosityGas, VolumeGas, ThermalConductivityGas, Chemical, Mixture


g = 9.81 # m / s2
R = 8.31446261815324 # m3 Pa / K mol
def K(C):
    return C+273.15
PI = 3.14159265358979323 # Why do we still remember this many digits

CASRN_water="7732-18-5"
CASRN_air="132259-10-0"

class GasMat:
    def __init__(self, code:str, M:float):
        '''
        code: CASRN code
        M: molar mass, in kg / mol
        '''

        self.code = code

        self.mu = ViscosityGas(CASRN=code)
        self.cp_base = HeatCapacityGas(CASRN=code, MW=M)
        # self.cp.method = "POLING_POLY"
        # self.chem = Chemical("air")
        self.k = ThermalConductivityGas(CASRN=code)

        self.M = M
        pass

    def vol(self, T_K, P):
        '''
        Assumes gas is ideal, and P is in Pa
        '''
        return R*T_K / P

    def beta(self, T_K):
        '''
        Assumes the gas is ideal
        '''
        return 1/T_K
    
    def cp(self, T):
        if self.code == CASRN_air:
            # Hacky solution because apparently there's no cp database with air in it by default
            # chem = Chemical(air, T=T)
            chem = Mixture(["N2", "O2"], zs=[0.787, 0.213], T=T)
            return chem.Cpg
        else:
            return self.cp_base(T)
    
    def alpha(self, T, P):
        # Thermal diffusivity
        return self.k(T, P)*self.vol(T, P)/self.cp(T)
    
    def rho(self, T, P):
        # Mass density
        return self.M / self.vol(T,P)
    
    def visc(self, T, P):
        # viscosity (m2/s)
        return self.mu(T, P)/self.rho(T, P)
    
    def Pr(self, T, P):
        # Prandtl number
        return self.visc(T, P)/self.alpha(T, P)
    
    def Ra(self, L, Ti, Ts, P):
        Ta = (Ti+Ts)/2
        return (g*self.beta(Ta)*abs(Ti-Ts)*L**3)/(self.visc(Ta, P)*self.alpha(Ta, P))
    

def Lc(r0, r1):
    '''
    Where r0 is the inner and r1 the outer radius of a concentric cylinder system.
    Note: Example 9.5's use of Eq. 9.59 appears to have a typo and has a value corresponding to the use of log10() not ln()
    '''
    log = np.log10 if USE_LOG10_TYPO else np.log
    return 2*((log(r1/r0))**(4/3)) / ((r0**-0.6 + r1**-0.6)**(5/3))

def Rac(gas:GasMat, r0, r1, T0, T1, P):
    T = (T0+T1)/2
    Lc_val = Lc(r0, r1)
    return (g*gas.beta(T)*abs(T0-T1)*(Lc_val**3)) / (gas.visc(T, P)*gas.alpha(T, P))

def keff(gas:GasMat, r0, r1, T0, T1, P):
    T = (T0+T1)/2
    Pr = gas.Pr(T, P)
    return 0.386*gas.k(T, P)*((Pr/(0.861+Pr))**0.25) * (Rac(gas, r0, r1, T0, T1, P)**0.25)
    

steam = GasMat(CASRN_water, 18.015e-3) # I stand corrected we didn't actually need to get this one
air = GasMat(CASRN_air, 28.96e-3)

# radii (m)
ri = 0.10 / 2
ro = 0.12 / 2
# rs = 0.14 / 2
rs = 0.15 / 2

USE_LOG10_TYPO = False # Set to false for authenticity, set to True to recreate the typo in Example 9.5 which calculates a value using log base 10 instead of the ln that it says it does.

Lca = Lc(ri, ro)
Lcb = Lc(ro, rs)
print(f"Lca: {Lca} m") # If USE_LOG10_TYPO=True, this is equal to Lc in 9.5
print(f"Lcb: {Lcb} m")


P = 101325 # Pa, 1 atm 

Ti = K(120) # inner steam
Ts = K(35) # outer whatever
To_INIT = (Ti+Ts)/2 # Start with average T
To = To_INIT

keffa = keff(air, ri, ro, Ti, To, P)
keffb = keff(air, ro, rs, To, Ts, P)
print(f"Starting keffa: {keffa} W/mK, keffb: {keffb}")
print(f"Starting To: {To:.3f}K")
for i in range(10):
    keffa = keff(air, ri, ro, Ti, To, P)
    keffb = keff(air, ro, rs, To, Ts, P)
    print(f"keffa: {keffa} W/mK, keffb: {keffb}")

    K_val = (keffa*np.log(rs/ro)) / (keffb*np.log(ro/ri))
    To = (Ts + K_val*Ti) / (1 + K_val)
    print(f"New To: {To:.3f}K")

qi_a = 2*PI*keffa*(Ti-To) / np.log(ro/ri)
qi_b = 2*PI*keffb*(To-Ts) / np.log(rs/ro)
print(f"Final q': {qi_a} W/m vs {qi_b}")







