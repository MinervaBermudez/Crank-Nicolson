# Heat Equation: Crank–Nicolson Finite-Difference Method

This repository contains a numerical solution of the **one-dimensional heat equation** using the Crank–Nicolson finite-difference method.

The project studies:

- The time evolution of a non-equilibrium initial temperature profile.
- The stationary solution associated with the boundary conditions.
- The convergence of the numerical solution towards thermal equilibrium.
- The distance to equilibrium using discrete $L^2$ and infinity norms.

The simulation was implemented in GNU Octave as part of a Numerical Analysis of Partial Differential Equations course.

---

## Mathematical model

The one-dimensional heat equation is

$$\frac{\partial u}{\partial t}=\frac{\partial^2u}{\partial x^2}.$$

Here:

- $u(x,t)$ represents the temperature at position $x$ and time $t$.
- $x\in[0,L]$ is the spatial coordinate.
- $t\in[0,T]$ is time.

The default domain and final time are

$$L=1,\qquad T=5.$$

The initial temperature profile is a Gaussian function centred at $x=0.5$:

$$u(x,0)=\exp\left[-20(x-0.5)^2\right].$$

---

## Boundary conditions

The implementation uses two different boundary conditions.

At the left boundary, a constant heat flux is imposed:

$$\frac{\partial u}{\partial x}(0,t)=-1.$$

At the right boundary, the temperature is fixed:

$$u(1,t)=0.$$

The left Neumann condition is discretised as

$$\frac{u_1-u_0}{\Delta x}=-1,$$

which gives

$$u_0=u_1+\Delta x.$$

This relation is used when reconstructing the solution at the left endpoint.

---

## Stationary solution

At equilibrium, the temperature no longer changes in time:

$$\frac{\partial u}{\partial t}=0.$$

Therefore, the equilibrium satisfies

$$\frac{d^2u_{\mathrm{eq}}}{dx^2}=0.$$

Its general form is

$$u_{\mathrm{eq}}(x)=Ax+B.$$

Applying the boundary conditions

$$u_{\mathrm{eq}}'(0)=-1,\qquad u_{\mathrm{eq}}(1)=0,$$

gives the unique equilibrium

$$u_{\mathrm{eq}}(x)=1-x.$$

The numerical solution is expected to approach this linear profile as time increases.

---

## Spatial and temporal discretisation

The spatial interval is divided into $N$ subintervals:

$$\Delta x=\frac{L}{N}.$$

The temporal interval is divided into $M$ steps:

$$\Delta t=\frac{T}{M}.$$

The default parameters are:

| Parameter | Value | Description |
| :--- | :--- | :--- |
| $L$ | $1$ | Length of the spatial domain. |
| $T$ | $5$ | Final simulation time. |
| $N$ | $50$ | Number of spatial subintervals. |
| $M$ | $100$ | Number of temporal steps. |
| $\Delta x$ | $0.02$ | Spatial step size. |
| $\Delta t$ | $0.05$ | Temporal step size. |

The second spatial derivative is approximated using centred finite differences:

$$\frac{\partial^2u}{\partial x^2}(x_i,t_n)\approx\frac{u_{i-1}^n-2u_i^n+u_{i+1}^n}{\Delta x^2}.$$

---

## Matrix formulation

The semi-discrete heat equation can be written as

$$\frac{d\mathbf{u}}{dt}=A\mathbf{u}+\mathbf{c},$$

where $\mathbf{u}$ contains the temperature at the interior grid points.

For the interior nodes, the finite-difference matrix has the form

$$A=\frac{1}{\Delta x^2}\begin{pmatrix}-1&1&0&\cdots&0\\1&-2&1&\cdots&0\\0&1&-2&\ddots&\vdots\\\vdots&\ddots&\ddots&\ddots&1\\0&\cdots&0&1&-2\end{pmatrix}.$$

The first diagonal coefficient is $-1/\Delta x^2$ instead of $-2/\Delta x^2$ because the Neumann boundary condition has been incorporated into the discretisation.

The constant boundary contribution is represented by

$$\mathbf{c}=\begin{pmatrix}1/\Delta x&0&\cdots&0\end{pmatrix}^{\mathsf T}.$$

---

## Crank–Nicolson method

The Crank–Nicolson method averages the spatial operator between two consecutive time levels:

$$\frac{\mathbf{u}^{n+1}-\mathbf{u}^n}{\Delta t}=\frac{1}{2}A\mathbf{u}^{n+1}+\frac{1}{2}A\mathbf{u}^n+\mathbf{c}.$$

Rearranging the equation gives

$$\left(I-\frac{\Delta t}{2}A\right)\mathbf{u}^{n+1}=\left(I+\frac{\Delta t}{2}A\right)\mathbf{u}^n+\Delta t\,\mathbf{c}.$$

At each time step, the program solves the linear system

```octave
b = right_matrix * u + g;
u = left_matrix \ b;
```

with

```octave
left_matrix  = I - (k/2)*A;
right_matrix = I + (k/2)*A;
```

The Crank–Nicolson method is second-order accurate in both space and time for sufficiently smooth solutions.

It is unconditionally stable for the linear heat equation, although the accuracy still depends on the selected values of $\Delta x$ and $\Delta t$.

---

## Time evolution

The initial Gaussian profile is not compatible with the stationary temperature distribution.

As the system evolves:

- Sharp spatial variations are smoothed by diffusion.
- The effect of the initial condition gradually disappears.
- The imposed boundary conditions determine the long-term temperature distribution.
- The solution approaches the equilibrium profile $u_{\mathrm{eq}}(x)=1-x$.

The numerical solution is represented both as an animated curve and as a three-dimensional surface over space and time.

---

## Distance to equilibrium

To quantify convergence, the numerical solution is compared with the equilibrium profile at every time step.

### Discrete $L^2$ distance

A mesh-weighted discrete approximation of the continuous $L^2$ norm is

$$d_2(t_n)=\left[\Delta x\sum_{i=1}^{N-1}\left(u_i^n-u_{\mathrm{eq}}(x_i)\right)^2\right]^{1/2}.$$

This quantity measures the global difference between the numerical solution and equilibrium.

In Octave, it can be calculated as

```octave
difference = u_save(:,j) - u_equilibrium(2:N)';
dist_l2(j) = sqrt(h * sum(difference.^2));
```

### Infinity distance

The discrete infinity norm is

$$d_\infty(t_n)=\max_i\left|u_i^n-u_{\mathrm{eq}}(x_i)\right|.$$

It measures the largest pointwise difference between the numerical solution and equilibrium.

In Octave:

```octave
difference = abs(u_save(:,j) - u_equilibrium(2:N)');
dist_inf(j) = max(difference);
```

Both distances are expected to decrease over time.

---

## Numerical results

The program generates the following visualisations:

1. Evolution of the temperature profile $u(x,t)$.
2. Three-dimensional surface of the solution over space and time.
3. Stationary equilibrium profile $u_{\mathrm{eq}}(x)=1-x$.
4. Evolution of the equilibrium profile under the numerical scheme.
5. Discrete $L^2$ distance to equilibrium.
6. Infinity-norm distance to equilibrium.

These plots illustrate the diffusive smoothing of the initial profile and the asymptotic convergence towards the unique stationary state.

---

## Requirements

- GNU Octave.
- A graphics toolkit supported by Octave.

No additional Octave packages are required.

---

## Running the simulation

Save the program as

```text
crank_nicolson_heat.m
```

and run it with

```bash
octave crank_nicolson_heat.m
```

The spatial and temporal resolution can be changed by modifying

```octave
N = 50;
M = 100;
```

The initial condition can be modified inside `u0_function`.
