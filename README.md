# Numerical Methods in MATLAB

Coursework from "Numerical Methods 1" and "Numerical Methods 2" at Warsaw University of Technology.

## Projects

### 1. Gauss-Seidel Solver (NM1)

Forward and backward Gauss-Seidel iterative methods for solving linear systems Ax = B. Includes convergence comparison, spectral radius analysis, and benchmarks up to 800×800 matrices.

```matlab
[X, iter, resHistory] = gaussSeidel(A, B, 1e-10, 1000);
```

### 2. Trapezoidal Integration (NM2)

2D composite trapezoidal rule for integration over a unit disk using polar coordinates. Tested on 10 functions with known exact solutions.

```matlab
I = approx_trapez_rule(@(x,y) exp(-(x.^2 + y.^2)), 32, 32);
```

### 3. Gill's RK4 ODE Solver (NM2)

Gill's variant of the 4th-order Runge-Kutta method for systems of 2 ODEs. Tested on 10 systems — 5 with exact solutions, 5 against ode45.

```matlab
[t, Y] = rk4_gill(@(t,Y) [Y(2); -4*Y(1)], 0, [1;0], 0.05, 100);
```

## Tech

MATLAB, numerical linear algebra, numerical integration, ODE solvers
