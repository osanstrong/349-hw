import json
import os.path

INPUT_PATH="input.json"
BASE_OUTPUT_PATH="output/output"
input = json.load(open(INPUT_PATH))

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
N_NODES = input["N_NODES"]
dr = R / (N_NODES - 0.5) # m, chosen so that final node (N-1) + 0.5 corresponds to the radius R, the exact surface  
T_0 = input["T_0_C"] # C
t_final = input["total_time_s"] # s

# Secondary setup (necessary dt)
SAFETY = input["SAFETY"] # Use smaller Fo than necessary
Bi = h*dr/k # Biot number, dimless
n_edge = N_NODES-0.5 # n of the last node, corresponds to radius of surface
n_inge = n_edge-0.5 # n of the border with second to last node (inside edge)
Fo = SAFETY * (1 / (3 * (n_inge**2/(n_edge**3 - n_inge**3)  +  Bi*n_edge**2/(n_edge**3 - n_inge**3))))
# print(f"Biot: {Bi}, Fo: {Fo}")
dt = Fo * dr**2 / alph
if "dt" in input and not input["dt"]=="auto":
    dt = input["dt"]
    Fo = dt*alph / (dr**2)
print(f"dt: {dt} s")

T_vals:list = [T_0 for i in range(N_NODES)] #
n_vals:list = [i+0.5 for i in range(N_NODES)] # The actual value n of each node (where radius = n*dr)


t_elapsed = 0
t_steps = []
T_vals_sets = [] #List of every T distribution at each timestep.
p = 0
num_timesteps = t_final / dt
num_prints = min(input["num_t_output"], num_timesteps)
print_freq = int(num_timesteps / num_prints)

def push():
    print(f"T({t_elapsed:2.3f}s) = {[f"{T:.3f}" for T in T_vals]}C")
    T_vals_sets.append(T_vals)
    t_steps.append(t_elapsed)

push()

while (t_elapsed < t_final):
    T_vals_next = [0 for i in range(N_NODES)] # Filler list
    
    if N_NODES > 1:
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
    else:
        T_nf = T_vals[0]
        T_vals_next[0] = T_nf + 3*Fo * (Bi*n_edge**2/(n_edge**3 - n_inge**3))*(T_INF-T_nf)

    
    t_elapsed += dt
    T_vals = T_vals_next
    p+=1
    
    if (p % print_freq == 0 or p < 10):
        push()
push()

output = {
    "input":input,
    "Biot (dr)":Bi,
    "Biot (D)":h*(2*R)/k,
    "t_steps":t_steps,
    "T_vals":T_vals_sets,
    "dt":dt
}

next_idx = 0
output_path = BASE_OUTPUT_PATH + ".json"
while os.path.isfile(output_path):
    output_path = BASE_OUTPUT_PATH + f"_{next_idx}.json"
    next_idx+=1
json.dump(output, open(output_path, "w"), indent=4)


