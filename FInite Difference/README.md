# Finite Difference Solver Project

This project implements the solution of the Poisson's equation in two-dimensional heterogeneous regions using the finite difference method.

## File Structure

- `FD.h`: This header file defines classes `FD_regular` and `FD_irregular`. Class `FD_regular` implements the solution of the Poisson's equation in the two-dimensional region [0,1]\*[0,1] using the finite difference method. Class `FD_irregular` implements the solution of the Poisson's equation in the two-dimensional region obtained by removing any circle from the region [0,1]\*[0,1] using the finite difference method. Both classes support three types of boundary conditions: Dirichlet, Neumann, and mixed boundary conditions.
- `test.cpp`: Test program. 
- `test.json`: Input parameter files for the test program.

## Test Programs

`test.cpp` is a test program that solves for three functions under different boundary conditions and regions. The functions are:

1. `u = e^(y + sin x)`
2. `u = e^(-xy)`
3. `u = sin(xy)`

The program changes the precision of the grid, allowing for analysis and comparison of the convergence rate at fixed points and the order of convergence of the error norm.

## Input Parameter Files

**test.json**

- `n`: The number of equal divisions of the unit square grid along the x and y directions. 
- `function_label`: The function number. Three functions have been designed for solution according to the problem requirements. Number 1 represents `u = e^(y + sin x)`, number 2 represents `u = e^(-xy)`, and number 3 represents `u = sin(xy)`. 
- `circle`: The information of the circle to be removed from the region. Three float numbers need to be inputted. The first represents the x-coordinate of the circle center, the second represents the y-coordinate of the circle center, and the third represents the radius of the circle. If the region is a unit square, set the radius to 0.
- `boundary_condition`: Input 1, 2, or 3 to represent the boundary condition. 1 represents Dirichlet, 2 represents Neumann, and 3 represents mixed.


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
