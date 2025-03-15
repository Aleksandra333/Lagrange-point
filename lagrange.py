import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
from scipy.optimize import fsolve


# Mass parameter: (Moon mass) / (Earth mass + Moon mass)
mu = 0.01215

def effective_potential(x, y):
    r1 = np.sqrt((x + mu)**2 + y**2)
    r2 = np.sqrt((x - 1 + mu)**2 + y**2)
    return 0.5*(x**2 + y**2) + (1 - mu) / r1 + mu / r2

def grad_U(X):
    x, y = X
    r1 = np.sqrt((x + mu)**2 + y**2)
    r2 = np.sqrt((x - 1 + mu)**2 + y**2)
    
    Ux = x - (1 - mu) * (x + mu) / r1**3 - mu * (x - 1 + mu) / r2**3
    Uy = y - (1 - mu) * y / r1**3 - mu * y / r2**3
    return np.array([Ux, Uy])

# The theoretical position:
t_x = 0.5 - mu
t_y = np.sqrt(3)/2

print(f"theoretical position:{t_x,t_y}")

t_L = np.array([t_x, t_y])
L, info, ier, mesg = fsolve(grad_U, t_L, full_output=True)

if ier != 1:
    print("fsolve did not converge:", mesg)
else:
    print("Lagrange point found at:", L)


# Equations in the rotating frame:
#   dx/dt = vx,
#   dy/dt = vy,
#   dvx/dt = U_x + 2 vy,
#   dvy/dt = U_y - 2 vx.
def deriv(t, state):
    x, y, vx, vy = state
    Ux, Uy = grad_U([x, y])
    
    # Include the Coriolis terms
    ax = Ux + 2 * vy
    ay = Uy - 2 * vx
    return [vx, vy, ax, ay]


perturbations = [
    np.array([0.001, 0.0]),
    np.array([0.0, 0.001]),
    np.array([0.0007, 0.0007]),
    np.array([-0.001, 0.001])
]

#store trajectories
trajectories = []

t_span = (0, 100)  # nondimensional time interval
t_eval = np.linspace(t_span[0], t_span[1], 5000)


for i, vel in enumerate(perturbations):
    # initial state: position at L, velocity = perturbation
    init_state = np.concatenate([L, vel])
    sol = solve_ivp(deriv, t_span, init_state, t_eval=t_eval, rtol=1e-9, atol=1e-12)
    trajectories.append(sol)


plt.figure(figsize=(8, 8))

# Earth at (-mu, 0)
#plt.plot(-mu, 0, 'bo', markersize=12, label="Earth")
# Moon at (1-mu, 0)
#plt.plot(1 - mu, 0, 'ko', markersize=8, label="Moon")
#Lagrange point
plt.plot(L[0], L[1], 'r*', markersize=15, label="Lagrange point (from fsolve)")

# Plot trajectories in the xy-plane.
for i, sol in enumerate(trajectories):
    plt.plot(sol.y[0], sol.y[1], label=f"Trajectory {i+1}")

plt.xlabel('x')
plt.ylabel('y')
plt.title('Trajectories near Lagrange point')
plt.legend()
plt.grid(True)
plt.axis('equal')
plt.show()


