# IVP Solver Project

This project implements eight different types of Initial Value Problem (IVP) solvers and uses them to solve the restricted three-body problem.

## The Restricted Three-Body Problem

The restricted three-body problem is a model for studying the motion of an earth-moon satellite or of an asteroid close enough to the earth to be strongly influenced by the earth and the sun. In this problem, we consider a system of three bodies: two heavy bodies, which move in circular orbits about their common center of mass, and a third body of negligible mass. 

In our formulation, we describe the position and velocity of the small body in a rotating coordinate system. The position coordinates of the small body are represented by the scalar variables `u1`, `u2`, and `u3`. Similarly, the velocity coordinates of the small body are represented by the scalar variables `u4`, `u5`, and `u6`. The masses of the two heavy bodies are represented by `1-mu` and `mu`.

The motion of the small body is governed by the following system of first-order differential equations:

```math
du1/dt = u4
du2/dt = u5
du3/dt = u6
du4/dt = 2*u5 + u1 - (1-mu)*(u1+mu)/r1^3 - mu*(u1-1+mu)/r2^3
du5/dt = -2*u4 + u2 - (1-mu)*u2/r1^3 - mu*u2/r2^3
du6/dt = -(1-mu)*u3/r1^3 - mu*u3/r2^3
```

where `r1 = sqrt((u1+mu)^2 + u2^2 + u3^2)` and `r2 = sqrt((u1-1+mu)^2 + u2^2 + u3^2)`.

In these equations, the gravitational influence of the two heavy bodies on the small body is evident. The small body is affected by both the centrifugal force (due to the rotating reference frame) and the gravitational forces from the two heavy bodies.

## Solvers Implemented

This project implements a variety of initial value problem (IVP) solvers:

- Adams-Bashforth methods of order 1, 2, 3, 4
- Adams-Moulton methods of order 2, 3, 4, 5
- Backward Differentiation Formulas of order 1, 2, 3, 4
- The classical Runge-Kutta (RK) method
- The Explicit Singly Diagonally Implicit Runge-Kutta (ESDIRK) method
- Gauss-Legendre RK methods of order 2, 3, 4, 5
- Fehlberg 4(5) embedded RK method
- Dormand-Prince 5(4) embedded RK method

Each of these solvers is designed to solve the restricted three-body problem, using different numerical methods to achieve the best possible accuracy and efficiency.

## File Structure

- `TI.h`: This is a header file that implements all of the IVP solvers mentioned in the Solvers Implemented section. Each call to an IVP solver will generate an object of (a derived class of) the `TimeIntegrator` class, instantiated by the singleton object `TimeIntegratorFactory`.
- `test.cpp`: Test program. 
- `test.json`: Input parameter files for the test program.

## Test Programs

The `test.cpp` program performs tests on orbits with two specific sets of initial values:

1. `(u1, u2, u3, u4, u5, u6) = (0.994, 0, 0, 0, -2.0015851063790825224, 0)`, which corresponds to a period `T1 = 17.06521656015796`.
2. `(u1, u2, u3, u4, u5, u6) = (0.879779227778, 0, 0, 0, -0.379677780949, 0)`, which corresponds to a period `T2 = 19.140540691377`.

For each orbit, the program runs orbital simulations and visualizations using all the IVP solvers implemented in the project. It then performs a sequence of grid refinement tests, where solution errors, convergence rates, and CPU time are reported. These tests demonstrate the convergence of each solver with a correct order of accuracy.

## Input Parameter Files

**test.json**

- `initial_value`: The initial value for the numerical algorithm. 
  - For the first orbit, the initial value should be set to `[0.994,0.0,0.0,0.0,-2.0015851063790825224,0.0]`
  - For the second orbit, the initial value should be set to `[0.879779227778,0.0,0.0,0.0,-0.379677780949,0.0]`

- `one_period_time`: The orbital period `T`. 
  - For the first orbit, the orbital period should be set to `17.06521656015796`
  - For the second orbit, the orbital period should be set to `19.140450691377`

- `one_period_step_number`: The number of numerical solution steps `n` per period, i.e., the number of iterations needed to iterate a period with the time-step size. Note that you can input a vector like `[10000,20000,40000,80000]` for calculating the convergence order. Please make sure that each element in the vector is twice the previous one. For the Fehlberg 4(5) embedded RK method and Dormand-Prince 5(4) embedded RK method, we will translate your input into the time-step size as the initial value due to their adaptive step size.

- `method`: The numerical solution algorithm. The eight numerical solution methods should be entered in top-to-bottom order as `Adams_Bashforth`, `Adams_Moulton`, `Backward_differentiation`, `Classical_RK`, `explicit_SDIRK`, `Gauss_Legendre_RK`, `Fehlberg_embedded_RK`, `Dormand_Prince_embedded_RK`.

- `order`: The order of accuracy `p` of the algorithm. The order of accuracy should be set for each method as follows:
  - `Adams_Bashforth` - 1,2,3,4
  - `Adams_Moulton` - 2,3,4,5
  - `Backward_differentiation` - 1,2,3,4
  - `Classical_RK` - 4
  - `explicit_SDIRK` - 4
  - `Gauss_Legendre_RK` - 2,4,6 (corresponding to s=1,2,3)
  - `Fehlberg_embedded_RK` - 4 (corresponding to 4(5))
  - `Dormand_Prince_embedded_RK` - 5 (corresponding to 5(4))

- `if_Richardson`: Whether to use Richardson extrapolation for order estimation. If "yes", the program will use Richardson extrapolation to give an estimate of the convergence order of the algorithm. If "no", Richardson extrapolation will not be used.

## Compilation and Execution

This project utilizes the `eigen3` library for solving systems of linear equations and uses the `jsoncpp` library for parameter input. These libraries are included in the project as `<eigen3/Eigen/...>` and `<jsoncpp/json/json.h>` respectively. Therefore, please ensure the correct file relationships while compiling.

The header file `TI.h` uses structure binding from C++17. Therefore, your computer needs to allow the `-std=c++17` used in `Makefile` for compilation.

You can compile the code by typing `make run` in the directory where the `Makefile` is located. This will create the executable file `test` respectively. Run to get the output for the test program.

```bash
# Compile
make run

# Run
./test
```

## Dependencies

This project depends on the following libraries:

1. [Eigen3]: A C++ template library for linear algebra. To install Eigen3, you can follow the instructions on their official website.

2. [JsonCpp]: A C++ library for interacting with JSON. You can install it using the following command:

    For Ubuntu:

    ```bash
    sudo apt-get install libjsoncpp-dev
    ```
