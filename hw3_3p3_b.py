import matplotlib.pyplot as plt
import numpy as np

# Don't change
T1inf = 40
h1 = 30

L = 0.004
k = 1.4


# Discrete change
h0_vals = [2, 65, 100]


# Continuous change
T0inf = np.linspace(-30, 0, 250)

def get_T1(giv_h0, giv_T0inf):
    return T1inf - (T1inf - giv_T0inf)/(h1 * (1/giv_h0 + L/k + 1/h1))

def get_T0(giv_h0, giv_T0inf):
    return giv_T0inf + (T1inf - giv_T0inf)/(giv_h0 * (1/giv_h0 + L/k + 1/h1))

for i in range(len(h0_vals)):
    h0 = h0_vals[i]
    c = f"C{i}"

    plt.plot(T0inf, get_T1(h0, T0inf), c=c, ls="--", label=f"Inner Surface, h0={h0} Wm⁻²K⁻1")
    plt.plot(T0inf, get_T0(h0, T0inf), c=c, ls="-", label=f"Outer Surface, h0={h0} Wm⁻²K⁻1")
    
plt.xlabel("Outside Ambient Air Temperature (⁰C)")
plt.ylabel("Temperature (⁰C)")
plt.legend()
plt.title("Problem 3.3b")
plt.show()