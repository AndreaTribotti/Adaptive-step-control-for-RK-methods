# Adaptive Step Control for Runge–Kutta Methods

This project studies **adaptive step-size control strategies for explicit Runge–Kutta (ERK) methods** applied to the numerical solution of **initial value problems (IVPs)** for ordinary differential equations.

The work was developed for the **Numerical Methods course** at the **University of Bologna**.

---

## Project Overview

The project is divided into two phases.

### Phase 1 – Adaptive Embedded Runge–Kutta Methods

We implement and test three embedded Runge–Kutta pairs with adaptive step control:

- **BS3(2)** – Bogacki–Shampine method  
- **ERK4(3)** – Embedded fourth/third-order Runge–Kutta  
- **DP5(4)** – Dormand–Prince method  

The methods are tested on:

- the **Brusselator system**
- two scalar ODE test problems.

#### How to run

To run a predefined experiment with fixed parameters:

```matlab
run_solver.m
```

This script runs the solvers and prints the results for a fixed configuration.

#### Interactive exploration

To explore the parameters interactively, open the MATLAB app:

```matlab
adaptive_steps_app.mlapp
```

The app allows the user to:

- select the ODE problem
- choose the Runge–Kutta method
- set absolute and relative tolerances
- modify the safety factor

It visualizes the solution, the step-size evolution and the errors behaviour.

---


### Phase 2 – Work–Precision Comparison

In the second phase, we compare **adaptive** and **fixed-step methods** on the **Brusselator system**.

**Adaptive methods:**
- BS3(2)
- ERK4(3)
- DP5(4)

**Fixed-step methods:**
- ERK3
- ERK4
- ERK4 with Richardson extrapolation

The performance comparison is based on a **work–precision diagram**, where **global error** is evaluated against the **number of function evaluations**.
To complement the global trends, a comparison table provides a **detailed performance snapshot for one representative configuration**.

#### Run the comparison

```matlab
run_full_comparison.m
```

This script generates the work–precision plot and the table used to compare the methods.

---

## Requirements

- MATLAB (with support for `.mlapp` apps)

---

## Detailed Report
For a full description of the methods, experiments, and results, see the full report: [*report_adaptive_step_control*](https://github.com/AndreaTribotti/Adaptive-step-control-for-RK-methods/blob/main/report_adaptive_step_control.pdf)  

---

## References

The adaptive step-size control strategy implemented in this project follows the framework described in:

E. Hairer, S. P. Nørsett, G. Wanner,  
[*Solving Ordinary Differential Equations I: Nonstiff Problems*,](https://github.com/AndreaTribotti/Adaptive-step-control-for-RK-methods/blob/main/Hairer_Norsett_Step_Control.pdf)  
Springer Series in Computational Mathematics.

---

## Authors

Andrea Tribotti  
Leonardo Cittadini  
Werther Solazzi  

