# Lagrange-point

## Lagrange Point Stability and Perturbed Trajectories

### Overview
This Python script investigates the stability of a Lagrange point in the circular restricted three-body problem (CRTBP) using numerical methods. It finds a Lagrange point by solving the effective potential gradient equations and then simulates perturbed trajectories near this point.

### Requirements
The script requires the following Python libraries:
- `numpy` for numerical computations
- `matplotlib` for plotting
- `scipy.integrate.solve_ivp` for numerical integration of differential equations
- `scipy.optimize.fsolve` for solving nonlinear equations

### Description
#### 1. **Mass Parameter**
The mass parameter `mu` is defined as the ratio of the Moon’s mass to the combined mass of the Earth and the Moon.

#### 2. **Effective Potential Calculation**
The script defines a function `effective_potential(x, y)` to compute the effective potential at a given position `(x, y)` in the rotating reference frame of the Earth-Moon system.

#### 3. **Gradient of Potential**
The function `grad_U(X)` computes the gradient of the effective potential, which is used to find equilibrium points.

#### 4. **Finding the Lagrange Point**
The function `fsolve()` from SciPy is used to solve `grad_U(X) = 0`, which gives the coordinates of a Lagrange point.

#### 5. **Equations of Motion**
The equations of motion in the rotating frame are defined in `deriv(t, state)`, incorporating the Coriolis effect:
- dx/dt = vx
- dy/dt = vy
- dvx/dt = Ux + 2vy
- dvy/dt = Uy - 2vx

#### 6. **Perturbed Orbits Simulation**
The script introduces small velocity perturbations near the Lagrange point and solves the equations using `solve_ivp()`. Several trajectories are computed with different initial velocity perturbations.

#### 7. **Plotting the Results**
The script generates a plot showing the computed Lagrange point and the perturbed trajectories.

### Expected Output
- The script prints the theoretically expected position of the Lagrange point and the numerically computed solution.
- A plot displaying:
  - The Earth and Moon’s positions (commented out but can be activated)
  - The computed Lagrange point
  - The perturbed trajectories in the rotating frame

### Example Output (Terminal)
```
theoretical position: (0.48785, 0.8660254037844386)
Lagrange point found at: [0.48785, 0.86603]
```

