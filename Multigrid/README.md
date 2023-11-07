# Multigrid Solver Project

This project implements the solution of one-dimensional and two-dimensional Poisson equations using various multigrid methods.

## File Structure

- `MC.h`: This header file provides the implementation of solving the Poisson equation using multigrid methods. Specifically, this file supports solving the Poisson equation for both **one-dimensional** cases within the interval [0,1] and **two-dimensional** cases within the area [0,1]\*[0,1]. Solutions support **Dirichlet**, **Neumann**, and **mixed boundary conditions**. The multigrid algorithm supports both **V-cycle** and **Full Multigrid (FMG)**. Restriction operators support both **full weighting** and **injection**. The interpolation operator provides support for both **linear** and **quadratic**.
- `test.cpp`: Test program. 
- `test.json`: Input parameter files for the test program.

## Test Programs

`test.cpp` is a test program. Specifically, it tests and compares the convergence rate of the residual norm at each iteration for both V-cycle and FMG in one-dimensional and two-dimensional scenarios. It also tests and compares the convergence rate of errors under different grid precisions for V-cycle and FMG. 

Moreover, for the one-dimensional scenario, the program tests the algorithm's ability to achieve the preset accuracy by gradually increasing the precision of the stopping criteria and determining the critical value at which the algorithm fails. 

For the two-dimensional scenario, the program calculates and analyzes the CPU time taken to solve the Poisson equation using V-cycle, FMG, and LU decomposition.

The test functions include:

**One-dimensional scenario:**

1. `y = e^(sin x)`
2. `y = sin x + 1`
3. `y = e^(-x)`

**Two-dimensional scenario:**

1. `y = e^(y + sin x)`
2. `y = sin(x + y) + 1`
3. `y = e^(-x + y)`

## Input Parameter Files

**test.json**

- `n`: The number of grid divisions. 
- `boundary_condition`: The boundary condition. You can choose between "Dirichlet", "Neumann", and "Mixed", which represent the corresponding boundary conditions.
- `restriction_operator`: The restriction operator. You can choose between "Full" and "Injection", which represent full-weighting and injection restriction operators.
- `interpolation_operator`: The interpolation operator. You can choose between "Linear" and "Quadratic", which represent first-order and second-order interpolation operators.
- `cycle`: The type of cycle. You can choose between "V" and "FMG", which represent V-cycle and Full Multigrid methods.
- `maximum_iteration_number`: The maximum number of iterations. The iteration stops when the number of iterations exceeds this value.
- `relative_accuracy`: The relative accuracy of the solution. The iteration stops when the relative accuracy of the solution is less than this value.
- `initial_guess`: The initial estimate of the solution. You can input "null vector" to represent a zero vector.

**test.cpp Parameter Settings**

In `test.cpp` file, you can specify certain parameters by modifying the following preprocessor definitions:

- `#define DIM`: By setting its value to 1 or 2, you can set the dimension parameter (`int dim`) of the template class to 1 or 2, respectively.
- `#define FUNCTION`: By setting its value to `Fun1`, `Fun2` or `Fun3`, you can specify the test function number as 1, 2 or 3, respectively. The corresponding expressions for these test functions can be found in the "Test Programs" section.


## Compilation and Execution

This project utilizes the `eigen3` library for solving systems of linear equations and uses the `jsoncpp` library for parameter input. These libraries are included in the project as `<eigen3/Eigen/...>` and `<jsoncpp/json/json.h>` respectively. Therefore, please ensure the correct file relationships while compiling.

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
