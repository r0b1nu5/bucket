import os
import numpy as np
from scipy.integrate import solve_ivp
#from scipy.integrate import odeint
import matplotlib.pyplot as plt
import xgi

def lorenz(t, xyz, sigma, rho, beta):
    x, y, z = xyz
    dxdt = sigma * (y - x)
    dydt = x * (rho - z) - y
    dzdt = x * y - beta * z
    return [dxdt, dydt, dzdt]

def coupled_lorenz(t, xyz, sigma, rho, beta, links, triangles, k2, k3):
    N = len(xyz) // 3  # Number of oscillators
    dxdt = np.zeros_like(xyz)
    
    for i in range(N):
        xi, yi, zi = xyz[i*3:i*3+3]
        dxdt[i*3:i*3+3] = lorenz(t, [xi, yi, zi], sigma, rho, beta)

    for i, j in links:
        xi, yi, zi = xyz[i*3:i*3+3]
        xj, yj, zj = xyz[j*3:j*3+3]
        # Coupling only directly influence the x variable
        dxdt[i*3:i*3+1] += k2*(xj-xi) # Coupled through the x variable
        dxdt[j*3:j*3+1] += k2*(xi-xj)

    for i, j, k in triangles:
        xi, yi, zi = xyz[i*3:i*3+3]
        xj, yj, zj = xyz[j*3:j*3+3]
        xk, yk, zk = xyz[k*3:k*3+3]
        # Coupling only directly influence the x variable
        dxdt[i*3:i*3+1] += k3*(xk*xj**2-xi**3) # Coupled through the x variable
        dxdt[j*3:j*3+1] += k3*(xk*xi**2-xj**3)
        dxdt[k*3:k*3+1] += k3*(xj*xi**2-xk**3)
    
    return dxdt

# Parameters
N = 5 # Number of oscillators
M = 10 # Number of different initial conditions
T = 3 # Simulation length for each initial condition
#time = np.linspace(0, T, 300) # Record data every 0.01 time unit
time = np.linspace(0,T,301)
sigma = 10.0
rho = 28.0
beta = 8.0 / 3.0
k2 = 1.0 # Pairwise coupling strength
k3 = 1.0 # Three-body coupling strength

# Random seed for reproducibility
np.random.seed(0)

# Generate random hypergraph (may want to switch to xgi.SimplicialComplex() for simplicial complexes)
H = xgi.random_hypergraph(N, [0.5, 0.5], seed=None)
H_int = xgi.convert_labels_to_integers(H, "label")
links = H_int.edges.filterby("size", 2).members()
triangles = H_int.edges.filterby("size", 3).members()
print(links)
print(triangles)

# Visualize the hypergraph
pos1 = xgi.random_layout(H)
pos2 = xgi.pairwise_spring_layout(H)
xgi.draw(H, pos2)

if os.path.exists('coupled_lorenz_solution.txt'):
    os.remove('coupled_lorenz_solution.txt')

f = open('coupled_lorenz_solution.txt','a')

for i in range(M):
    # Initial conditions (may want to control the box size)
    xyz0 = 1e-1*np.random.randn(N*3) + 20*np.ones(N*3)
    #sol = odeint(coupled_lorenz, xyz0, time, args=(sigma, rho, beta, links, triangles, k2, k3))
    sol = solve_ivp(coupled_lorenz, [0, T], xyz0, t_eval=time, args=(sigma, rho, beta, links, triangles, k2, k3), max_step=1e-3)# rtol=1e-9, atol=1e-9)
    print(np.shape(sol.y.T))
    np.savetxt(f, sol.y.T)

f.close()

# Plot the last solution
t = sol.t
xyz = sol.y

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

for i in range(N):
    ax.plot(xyz[i*3], xyz[i*3+1], xyz[i*3+2])

ax.set_xlabel('X')
ax.set_ylabel('Y')
ax.set_zlabel('Z')
plt.show()
