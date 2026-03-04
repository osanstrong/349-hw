
# Copper Material Properties
rho = 8933 # kg / m^3
c_p = 385 # J / kg*K
k = 401 # W / m*K
alph = k / (rho*c_p) # (m^2 / s)

# Fluid properties
h = 49.96 # W/m^2*K
T_INF = 25 # C

# Experiment Setup
R = 0.01 # m
N_NODES = 40
dr = R / (N_NODES - 0.5) # m, chosen so that final node (N-1) + 0.5 corresponds to the radius R, the exact surface  
T_0 = 84 # C
t_final = 95 # s

# Secondary setup (necessary dt)
SAFETY = 0.75 # Use smaller Fo than necessary
Bi = h*dr/k # Biot number, dimless
n_edge = N_NODES-0.5 # n of the last node, corresponds to radius of surface
n_inge = n_edge-0.5 # n of the border with second to last node (inside edge)
Fo = SAFETY * (1 / (3 * (n_inge**2/(n_edge**3 - n_inge**3)  +  Bi*n_edge**2/(n_edge**3 - n_inge**3))))
# print(f"Biot: {Bi}, Fo: {Fo}")
dt = Fo * dr**2 / alph
# print(f"dt: {dt} s")
print(f"Fourier: {Fo} vs {alph * dt / (dr**2)}")

T_vals:list = [T_0 for i in range(N_NODES)] #
n_vals:list = [i+0.5 for i in range(N_NODES)] # The actual value n of each node (where radius = n*dr)


t_elapsed = 0
T_vals_sets = [] #List of every T distribution at each timestep.
p = 0
num_timesteps = t_final / dt
num_prints = 15
print_freq = int(num_timesteps / num_prints)
while (t_elapsed < t_final):
    T_vals_next = [0 for i in range(N_NODES)] # Filler list
    
    # Center node
    T_n0 = T_vals[0]
    T_n1 = T_vals[1]
    n0 = n_vals[0]
    nh = n0 + 0.5
    T_vals_next[0] = T_n0 + Fo * (nh**2 / n0**2) * (T_n1-T_n0) 

    for i in range(1, N_NODES-1): # Interior nodes, not special cases
        n = n_vals[i]
        T_n = T_vals[i]
        T_i = T_vals[i-1]
        T_o = T_vals[i+1]

        T_vals_next[i] = T_n + Fo * ((n+0.5)**2*(T_o-T_n) - (n-0.5)**2*(T_n-T_i)) / n**2
        # print(f"new T_{i}: {T_vals_next[i]}")

    # Exterior node
    T_nf = T_vals[-1]
    T_vals_next[-1] = T_nf + 3*Fo * ((n_inge**2/(n_edge**3 - n_inge**3))*(T_vals[-2]-T_nf) + (Bi*n_edge**2/(n_edge**3 - n_inge**3))*(T_INF-T_nf))
    

    T_vals_sets.append(T_vals)
    t_elapsed += dt
    T_vals = T_vals_next
    p+=1
    
    if (p % print_freq == 0 or t_elapsed >= t_final): print(f"T({t_elapsed:2.3f}s) = {[f"{T:.3f}" for T in T_vals_next]}C")


