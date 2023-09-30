#include <iostream>
#include <cmath>
#include <algorithm>
#include <vector>
#include <chrono>
#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/Sparse>
#include <eigen3/Eigen/Core>


template <int dim>
class Function
{
public:
    /*
    Univariate Function
    */
    virtual double operator()(double _x) = 0;       // Return the value of function u at _x
    virtual double diff(double _x) { return 0; }    // Return the first order derivative of function u at _x
    virtual double Laplace(double _x) { return 0; } // Return the value of function u at _x after the -Laplace operator is applied

    /*
    Bivariate Function
    */
    virtual double operator()(double _x, double _y) = 0;       // Return the value of function u at (_x,_y)
    virtual double diff_x(double _x, double _y) { return 0; }  // Return the partial derivative of function u with respect to x at (_x,_y)
    virtual double diff_y(double _x, double _y) { return 0; }  // Return the partial derivative of function u with respect to y at (_x,_y)
    virtual double Laplace(double _x, double _y) { return 0; } // Return the value of function u at (_x,_y) after the -Laplace operator is applied
};

/*
    Weighted Jacobi Iteration
*/
template <int dim>
class Weighted_Jacobi
{
private:
    Eigen::SparseMatrix<double> A; 
    Eigen::VectorXd b;
    Eigen::VectorXd init_x; // Initial value for iteration
    int iter_Num;           // Number of iterations
    double w;               // Iteration weight
public:
    Weighted_Jacobi(Eigen::SparseMatrix<double> _A, Eigen::VectorXd _b, Eigen::VectorXd _init_x, int _iter_Num, double _w) : A(_A), b(_b), init_x(_init_x), iter_Num(_iter_Num), w(_w) {}

    Eigen::VectorXd solve()
    {

        Eigen::SparseMatrix<double> D(A.rows(), A.cols());
        Eigen::SparseMatrix<double> D_inverse(A.rows(), A.cols());

        for (int i = 0; i < A.rows(); i++)
        {
            D.insert(i, i) = A.coeff(i, i);
            D_inverse.insert(i, i) = 1.0 / A.coeff(i, i);
        }

        Eigen::SparseMatrix<double> T = D_inverse * (D - A);
        Eigen::VectorXd c = D_inverse * b;
        Eigen::VectorXd x = init_x;

        for (int k = 0; k < iter_Num; k++)
        {
            Eigen::VectorXd x_star = T * x + c;
            x = (1.0 - w) * x + w * x_star;
        }

        return x;
    }
};

/*
    interpolation and estriction
*/
template <int dim>
class Transformation
{
private:
    Eigen::VectorXd v; // Vector to be processed
    double h;          // Grid width
public:
    Transformation(Eigen::VectorXd _v, double _h) : v(_v), h(_h) {}

    // interpolation operator
    Eigen::VectorXd Interpolation(std::string inter, std::string boundary)
    {
        if (dim == 1)
        {
            if (inter == "Linear")
            {
                int n = static_cast<int>(2.0 / h + 1.0);
                Eigen::SparseMatrix<double> I(n, (n + 1) / 2);

                for (int i = 1; i < (n - 1) / 2; i++)
                {
                    I.insert(2 * i, i) = 1.0;
                    I.insert(2 * i - 1, i) = 0.5;
                    I.insert(2 * i + 1, i) = 0.5;
                }
                I.insert(0, 0) = 1.0;
                I.insert(1, 0) = 0.5;
                I.insert(n-1, (n-1)/2) = 1.0;
                I.insert(n-2, (n-1)/2) = 0.5;

                return I * v;
            }
            if (inter == "Quadratic")
            {
                int n = static_cast<int>(2.0 / h + 1.0);
                Eigen::SparseMatrix<double> I(n, (n + 1) / 2);

                for (int i = 0; i < (n + 1) / 2; i++)
                {
                    I.insert(2 * i, i) = 1.0;
                }
                for (int i = 1; i < (n - 3) / 2; i++)
                {
                    I.insert(2 * i + 1, i) = 9.0/16.0;
                    I.insert(2 * i + 1, i - 1) = -1.0/16.0;
                    I.insert(2 * i + 1, i + 1) = 9.0/16.0;
                    I.insert(2 * i + 1, i + 2) = -1.0/16.0;
                }
                I.insert(1, 1) = 0.5;
                I.insert(1, 0) = 0.5;
                I.insert(n-2, (n-3)/2) = 0.5;
                I.insert(n-2, (n-1)/2) = 0.5;

                return I * v;
            }
        }
        if (dim == 2)
        {
            if (inter == "Linear")
            {
                int n_old = static_cast<int>(1.0 / h + 1.0);
                int n_new = static_cast<int>(2.0 / h + 1.0);
                Eigen::VectorXd v_new(n_new*n_new);

                for (int i_new = 0; i_new < n_new*n_new; i_new++)
                {
                    int Label_x_new = i_new % n_new;
                    int Label_y_new = i_new / n_new;
                    if (Label_x_new % 2 == 0 && Label_y_new % 2 == 0) 
                    {
                        v_new(i_new) = v(Label_x_new/2 + Label_y_new/2*n_old);
                    }
                    if (Label_x_new % 2 == 1 && Label_y_new % 2 == 0) 
                    {
                        v_new(i_new) = 0.5*v((Label_x_new+1)/2 + Label_y_new/2*n_old) + 0.5*v((Label_x_new-1)/2 + Label_y_new/2*n_old);
                    }
                    if (Label_x_new % 2 == 0 && Label_y_new % 2 == 1) 
                    {
                        v_new(i_new) = 0.5*v(Label_x_new/2 + (Label_y_new+1)/2*n_old) + 0.5*v(Label_x_new/2 + (Label_y_new-1)/2*n_old);
                    }
                    if (Label_x_new % 2 == 1 && Label_y_new % 2 == 1) 
                    {
                        v_new(i_new) = 0.25*v((Label_x_new+1)/2 + (Label_y_new+1)/2*n_old) + 0.25*v((Label_x_new+1)/2 + (Label_y_new-1)/2*n_old)
                                     + 0.25*v((Label_x_new-1)/2 + (Label_y_new+1)/2*n_old) + 0.25*v((Label_x_new-1)/2 + (Label_y_new-1)/2*n_old);
                    }
                }
               
                return v_new;
            }
            if (inter == "Quadratic")
            {
                int n_old = static_cast<int>(1.0 / h + 1.0);
                int n_new = static_cast<int>(2.0 / h + 1.0);
                Eigen::VectorXd v_new(n_new*n_new);

                for (int i_new = 0; i_new < n_new*n_new; i_new++)
                {
                    int Label_x_new = i_new % n_new;
                    int Label_y_new = i_new / n_new;
                    if (Label_x_new % 2 == 0 && Label_y_new % 2 == 0) 
                    {
                        v_new(i_new) = v(Label_x_new/2 + Label_y_new/2*n_old);
                    }
                    if (Label_x_new % 2 == 1 && Label_y_new % 2 == 0) 
                    {
                        if (Label_x_new == 1 || Label_x_new == n_new - 2)
                        {
                            v_new(i_new) = 0.5*v((Label_x_new+1)/2 + Label_y_new/2*n_old) + 0.5*v((Label_x_new-1)/2 + Label_y_new/2*n_old);
                        }
                        else
                        {
                            v_new(i_new) = 9.0/16.0*v((Label_x_new+1)/2 + Label_y_new/2*n_old) + 9.0/16.0*v((Label_x_new-1)/2 + Label_y_new/2*n_old)
                                           -1.0/16.0*v((Label_x_new+1)/2 + Label_y_new/2*n_old + 1) - 1.0/16.0*v((Label_x_new-1)/2 + Label_y_new/2*n_old - 1);
                        }
                    }
                    if (Label_x_new % 2 == 0 && Label_y_new % 2 == 1) 
                    {
                        if (Label_y_new == 1 || Label_y_new == n_new - 2)
                        {
                            v_new(i_new) = 0.5*v(Label_x_new/2 + (Label_y_new+1)/2*n_old) + 0.5*v(Label_x_new/2 + (Label_y_new-1)/2*n_old);
                        }
                        else
                        {
                            v_new(i_new) = 9.0/16.0*v(Label_x_new/2 + (Label_y_new+1)/2*n_old) + 9.0/16.0*v(Label_x_new/2 + (Label_y_new-1)/2*n_old);
                                           -1.0/16.0*v(Label_x_new/2 + (Label_y_new+1)/2*n_old + n_old) - 1.0/16.0*v(Label_x_new/2 + (Label_y_new-1)/2*n_old - n_old);
                        }
                    }
                    if (Label_x_new % 2 == 1 && Label_y_new % 2 == 1) 
                    {
                        if (Label_y_new == 1 || Label_y_new == n_new - 2 || Label_x_new == 1 || Label_x_new == n_new - 2)
                        {
                            v_new(i_new) = 0.25*v((Label_x_new+1)/2 + (Label_y_new+1)/2*n_old) + 0.25*v((Label_x_new+1)/2 + (Label_y_new-1)/2*n_old)
                                         + 0.25*v((Label_x_new-1)/2 + (Label_y_new+1)/2*n_old) + 0.25*v((Label_x_new-1)/2 + (Label_y_new-1)/2*n_old);
                        }
                        else
                        {
                            v_new(i_new) = 9.0/32.0*v((Label_x_new+1)/2 + (Label_y_new+1)/2*n_old) + 9.0/32.0*v((Label_x_new+1)/2 + (Label_y_new-1)/2*n_old)
                                         + 9.0/32.0*v((Label_x_new-1)/2 + (Label_y_new+1)/2*n_old) + 9.0/32.0*v((Label_x_new-1)/2 + (Label_y_new-1)/2*n_old)
                                         - 1.0/64.0*v((Label_x_new+1)/2 + (Label_y_new+1)/2*n_old + 1) - 1.0/64.0*v((Label_x_new+1)/2 + (Label_y_new+1)/2*n_old + n_old)
                                         - 1.0/64.0*v((Label_x_new+1)/2 + (Label_y_new-1)/2*n_old + 1) - 1.0/64.0*v((Label_x_new+1)/2 + (Label_y_new-1)/2*n_old - n_old)
                                         - 1.0/64.0*v((Label_x_new-1)/2 + (Label_y_new+1)/2*n_old - 1) - 1.0/64.0*v((Label_x_new-1)/2 + (Label_y_new+1)/2*n_old + n_old)
                                         - 1.0/64.0*v((Label_x_new-1)/2 + (Label_y_new-1)/2*n_old - 1) - 1.0/64.0*v((Label_x_new-1)/2 + (Label_y_new-1)/2*n_old - n_old);
                        }
                    }
                }
               
                return v_new;
            }
        }
        Eigen::VectorXd vec;
        return vec;
    }

    // restriction operator
    Eigen::VectorXd Restriction(std::string res, std::string boundary)
    {
        if (dim == 1)
        {
            if (res == "Full")
            {
                int n = static_cast<int>(1.0 / (2.0 * h) + 1.0);

                Eigen::SparseMatrix<double> I(n, n * 2 - 1);

                for (int i = 1; i < n - 1; i++)
                {
                    I.insert(i, 2 * i) = 0.5;
                    I.insert(i, 2 * i - 1) = 0.25;
                    I.insert(i, 2 * i + 1) = 0.25;
                }
                I.insert(0, 0) = 1.0;
                I.insert(n - 1, n * 2 - 2) = 1.0;           

                return I * v;
            }
            if (res == "Injection")
            {
                int n = static_cast<int>(1.0 / (2.0 * h) + 1.0);

                Eigen::SparseMatrix<double> I(n, n * 2 - 1);

                for (int i = 0; i < n ; i++)
                {
                    I.insert(i, 2 * i) = 1.0;
                }

                return I * v;
            }
        }
        if (dim == 2)
        {
            if (res == "Full")
            {
                int n_old = static_cast<int>(1.0 / h + 1.0);
                int n_new = static_cast<int>(1.0 / (2.0 * h) + 1.0);
                Eigen::VectorXd v_new(n_new*n_new);

                for (int i_new = 0; i_new < n_new*n_new; i_new++)
                {
                    int Label_x_new = i_new % n_new;
                    int Label_y_new = i_new / n_new;
                    if (i_new % n_new == 0 || i_new % n_new == n_new-1 || i_new / n_new == 0 || i_new / n_new == n_new-1)
                    {
                        v_new(i_new) = v(Label_x_new*2+Label_y_new*2*n_old);
                    }
                    else
                    {
                        v_new(i_new) = 1/4.0*v(Label_x_new*2 + Label_y_new*2*n_old)
                                     + 1/8.0*v(Label_x_new*2 + Label_y_new*2*n_old + 1) + 1/8.0*v(Label_x_new*2 + Label_y_new*2*n_old - 1)
                                     + 1/8.0*v(Label_x_new*2 + Label_y_new*2*n_old + n_old) + 1/8.0*v(Label_x_new*2 + Label_y_new*2*n_old - n_old)
                                     + 1/16.0*v(Label_x_new*2 + Label_y_new*2*n_old + n_old + 1) + 1/8.0*v(Label_x_new*2 + Label_y_new*2*n_old + n_old - 1)
                                     + 1/16.0*v(Label_x_new*2 + Label_y_new*2*n_old - n_old + 1) + 1/8.0*v(Label_x_new*2 + Label_y_new*2*n_old - n_old - 1);
                    }                
                }

                return v_new;
            }
            if (res == "Injection")
            {
                int n_old = static_cast<int>(1.0 / h + 1.0);
                int n_new = static_cast<int>(1.0 / (2.0 * h) + 1.0);
                Eigen::VectorXd v_new(n_new*n_new);

                for (int i_new = 0; i_new < n_new*n_new; i_new++)
                {
                    int Label_x_new = i_new % n_new;
                    int Label_y_new = i_new / n_new;
                    v_new(i_new) = v(Label_x_new*2+Label_y_new*2*n_old);        
                }

                return v_new;
            }
        }
        Eigen::VectorXd vec;
        return vec;
    }
};

/*
    Construct the coefficient matrix A_h and the initial right-hand side vector init_f
*/
template <int dim>
class Coefficient
{
private:
    Function<dim> &u;     // Given function
    double h;             // Grid width
    std::string boundary; // Boundary conditions
public:
    Coefficient(Function<dim> &_u, std::string _boundary, double _h) : u(_u), boundary(_boundary), h(_h) {}

    // Construct matrix A_h
    Eigen::SparseMatrix<double> get_Ah()
    {
        if (dim == 1)
        {
            if (boundary == "Dirichlet")
            {
                int n = static_cast<int>(1.0 / h + 1.0);
                Eigen::SparseMatrix<double> Ah(n, n);
                std::vector<Eigen::Triplet<double> > tripletlist;

                for (int i = 1; i < n - 1; ++i)
                {
                    tripletlist.push_back(Eigen::Triplet<double>(i, i, 2.0/(h*h)));
                    tripletlist.push_back(Eigen::Triplet<double>(i, i-1, -1.0/(h*h)));
                    tripletlist.push_back(Eigen::Triplet<double>(i, i+1, -1.0/(h*h)));
                }
                tripletlist.push_back(Eigen::Triplet<double>(0, 0, 1.0));
                tripletlist.push_back(Eigen::Triplet<double>(n-1, n-1, 1.0));
                Ah.setFromTriplets(tripletlist.begin(), tripletlist.end());
                Ah.makeCompressed();

                return Ah;
            }
            if (boundary == "Neumann")
            {
                int n = static_cast<int>(1.0 / h + 1.0);
                Eigen::SparseMatrix<double> Ah(n, n);
                std::vector<Eigen::Triplet<double> > tripletlist;

                for (int i = 1; i < n - 1; ++i)
                {
                    tripletlist.push_back(Eigen::Triplet<double>(i, i, 2.0/(h*h)));
                    tripletlist.push_back(Eigen::Triplet<double>(i, i-1, -1.0/(h*h)));
                    tripletlist.push_back(Eigen::Triplet<double>(i, i+1, -1.0/(h*h)));
                }
                tripletlist.push_back(Eigen::Triplet<double>(0, 0, 1.0/h));
                tripletlist.push_back(Eigen::Triplet<double>(0, 1, 2.0/(h)));
                tripletlist.push_back(Eigen::Triplet<double>(0, 2, -1.0/(2.0*h)));
                tripletlist.push_back(Eigen::Triplet<double>(n-1, n-1, 3.0/(2.0*h)));
                tripletlist.push_back(Eigen::Triplet<double>(n-1, n-2, -2.0/(h)));
                tripletlist.push_back(Eigen::Triplet<double>(n-1, n-3, 1.0/(2.0*h)));
                Ah.setFromTriplets(tripletlist.begin(), tripletlist.end());
                Ah.makeCompressed();

                return Ah;
            }
            if (boundary == "Mixed")
            {
                int n = static_cast<int>(1.0 / h + 1.0);
                Eigen::SparseMatrix<double> Ah(n, n);
                std::vector<Eigen::Triplet<double> > tripletlist;

                for (int i = 1; i < n - 1; ++i)
                {
                    tripletlist.push_back(Eigen::Triplet<double>(i, i, 2.0/(h*h)));
                    tripletlist.push_back(Eigen::Triplet<double>(i, i-1, -1.0/(h*h)));
                    tripletlist.push_back(Eigen::Triplet<double>(i, i+1, -1.0/(h*h)));
                }
                tripletlist.push_back(Eigen::Triplet<double>(0, 0, 1.0));
                tripletlist.push_back(Eigen::Triplet<double>(n-1, n-1, 3.0/(2.0*h)));
                tripletlist.push_back(Eigen::Triplet<double>(n-1, n-2, -2.0/(h)));
                tripletlist.push_back(Eigen::Triplet<double>(n-1, n-3, 1.0/(2.0*h)));
                Ah.setFromTriplets(tripletlist.begin(), tripletlist.end());
                Ah.makeCompressed();

                return Ah;
            }
        }
        if (dim == 2)
        {
            if (boundary == "Dirichlet")
            {
                int n = static_cast<int>(1.0 / h + 1.0);
                Eigen::SparseMatrix<double> Ah(n*n, n*n);
                std::vector<Eigen::Triplet<double> > tripletlist;

                for (int i = 0; i < n*n; i++)
		        {
                    if (i % n == 0 || i % n == n-1 || i / n == 0 || i / n == n-1)
                    {
                        tripletlist.push_back(Eigen::Triplet<double>(i, i, 1.0));
                    }
                    else
                    {
                        tripletlist.push_back(Eigen::Triplet<double>(i, i, 4.0/(h*h)));
                        tripletlist.push_back(Eigen::Triplet<double>(i, i-1, -1.0/(h*h)));
                        tripletlist.push_back(Eigen::Triplet<double>(i, i+1, -1.0/(h*h)));
                        tripletlist.push_back(Eigen::Triplet<double>(i, i-n, -1.0/(h*h)));
                        tripletlist.push_back(Eigen::Triplet<double>(i, i+n, -1.0/(h*h)));
                    }
                }
                Ah.setFromTriplets(tripletlist.begin(), tripletlist.end());
                Ah.makeCompressed();

                return Ah;
            }
            if (boundary == "Neumann")
            {
                int n = static_cast<int>(1.0 / h + 1.0);
                Eigen::SparseMatrix<double> Ah(n*n, n*n);
                std::vector<Eigen::Triplet<double> > tripletlist;

                for (int i = 0; i < n*n; i++)
		        {
                    if (i == 0 || i == n-1 || i == n*n-n || i == n*n-1)
                    {
                        tripletlist.push_back(Eigen::Triplet<double>(i, i, 1.0/(h*h)));
                    }
                    else
                    {
                        if (i % n == 0 || i % n == n-1 || i / n == 0 || i / n == n-1)
                        {
                            if (i % n == 0) // Left Boundary
                            {
                                tripletlist.push_back(Eigen::Triplet<double>(i, i, -3.0/(2.0*h*h)));
                                tripletlist.push_back(Eigen::Triplet<double>(i, i+1, 2.0/(h*h)));
                                tripletlist.push_back(Eigen::Triplet<double>(i, i+2, -1.0/(2.0*h*h)));
                            }
                            if (i % n == n-1) // Right Boundary
                            {
                                tripletlist.push_back(Eigen::Triplet<double>(i, i, 3.0/(2.0*h*h)));
                                tripletlist.push_back(Eigen::Triplet<double>(i, i-1, -2.0/(h*h)));
                                tripletlist.push_back(Eigen::Triplet<double>(i, i-2, 1.0/(2.0*h*h)));
                            }
                            if (i / n == 0) // Lower Boundary
                            {
                                if (i == 1)
                                {
                                    tripletlist.push_back(Eigen::Triplet<double>(i, i, 1.0/(h*h)));
                                }
                                else
                                {
                                    tripletlist.push_back(Eigen::Triplet<double>(i, i, -3.0/(2.0*h*h)));
                                    tripletlist.push_back(Eigen::Triplet<double>(i, i+n, 2.0/(h*h)));
                                    tripletlist.push_back(Eigen::Triplet<double>(i, i+2*n, -1.0/(2.0*h*h)));
                                } 
                            }
                            if (i / n == n-1) // Upper Boundary
                            {
                                tripletlist.push_back(Eigen::Triplet<double>(i, i, 3.0/(2.0*h*h)));
                                tripletlist.push_back(Eigen::Triplet<double>(i, i-n, -2.0/(h*h)));
                                tripletlist.push_back(Eigen::Triplet<double>(i, i-2*n, 1.0/(2.0*h*h)));
                            }
                        }
                        else
                        {
                            tripletlist.push_back(Eigen::Triplet<double>(i, i, 4.0/(h*h)));
                            tripletlist.push_back(Eigen::Triplet<double>(i, i-1, -1.0/(h*h)));
                            tripletlist.push_back(Eigen::Triplet<double>(i, i+1, -1.0/(h*h)));
                            tripletlist.push_back(Eigen::Triplet<double>(i, i-n, -1.0/(h*h)));
                            tripletlist.push_back(Eigen::Triplet<double>(i, i+n, -1.0/(h*h)));
                        }
                    }
		        }
                Ah.setFromTriplets(tripletlist.begin(), tripletlist.end());
                Ah.makeCompressed();

                return Ah;
            }
            if (boundary == "Mixed")
            {
                int n = static_cast<int>(1.0 / h + 1.0);
                Eigen::SparseMatrix<double> Ah(n*n, n*n);
                std::vector<Eigen::Triplet<double> > tripletlist;

                for (int i = 0; i < n*n; i++)
		        {
                    if (i == 0 || i == n-1 || i == n*n-n || i == n*n-1)
                    {
                        tripletlist.push_back(Eigen::Triplet<double>(i, i, 1.0));
                    }
                    else
                    {
                        if (i % n == 0 || i % n == n-1 || i / n == 0 || i / n == n-1)
                        {
                            if (i % n == 0) // Left Boundary
                            {
                                tripletlist.push_back(Eigen::Triplet<double>(i, i, -3.0/(2.0*h)));
                                tripletlist.push_back(Eigen::Triplet<double>(i, i+1, 2.0/(h)));
                                tripletlist.push_back(Eigen::Triplet<double>(i, i+2, -1.0/(2.0*h)));
                            }
                            if (i % n == n-1) // Right Boundary
                            {
                                tripletlist.push_back(Eigen::Triplet<double>(i, i, 3.0/(2.0*h)));
                                tripletlist.push_back(Eigen::Triplet<double>(i, i-1, -2.0/(h)));
                                tripletlist.push_back(Eigen::Triplet<double>(i, i-2, 1.0/(2.0*h)));
                            }
                            if (i / n == 0) // Lower Boundary
                            {
                                tripletlist.push_back(Eigen::Triplet<double>(i, i, 1.0));
                            }
                            if (i / n == n-1) // Upper Boundary
                            {
                                tripletlist.push_back(Eigen::Triplet<double>(i, i, 1.0));
                            }
                        }
                        else
                        {
                            tripletlist.push_back(Eigen::Triplet<double>(i, i, 4.0/(h*h)));
                            tripletlist.push_back(Eigen::Triplet<double>(i, i-1, -1.0/(h*h)));
                            tripletlist.push_back(Eigen::Triplet<double>(i, i+1, -1.0/(h*h)));
                            tripletlist.push_back(Eigen::Triplet<double>(i, i-n, -1.0/(h*h)));
                            tripletlist.push_back(Eigen::Triplet<double>(i, i+n, -1.0/(h*h)));
                        }
                    }
		        }
                Ah.setFromTriplets(tripletlist.begin(), tripletlist.end());
                Ah.makeCompressed();

                return Ah;
            }
        }
        Eigen::SparseMatrix<double> vec;
        return vec;
    }

    Eigen::VectorXd get_initf()
    {
        if (dim == 1)
        {
            if (boundary == "Dirichlet")
            {
                int n = static_cast<int>(1.0 / h + 1.0);
                Eigen::VectorXd f_init(n);

                f_init(0) = u(0);
                f_init(n - 1) = u(1);

                for (int i = 1; i < n - 1; i++)
                {
                    f_init(i) = u.Laplace(i * h);
                }

                return f_init;
            }
            if (boundary == "Neumann")
            {
                int n = static_cast<int>(1.0 / h + 1.0);
                Eigen::VectorXd f_init(n);

                f_init(0) = u.diff(0) + (3.0*u(0))/(2.0*h) + (1.0*u(0))/h;
                f_init(n - 1) = u.diff(1);

                for (int i = 1; i < n - 1; i++)
                {
                    f_init(i) = u.Laplace(i * h);
                }

                return f_init;
            }
            if (boundary == "Mixed")
            {
                int n = static_cast<int>(1.0 / h + 1.0);
                Eigen::VectorXd f_init(n);

                f_init(0) = u(0);
                f_init(n - 1) = u.diff(1);

                for (int i = 1; i < n - 1; i++)
                {
                    f_init(i) = u.Laplace(i * h);
                }

                return f_init;
            }
        }
        if (dim == 2)
        {
            if (boundary == "Dirichlet")
            {
                int n = static_cast<int>(1.0 / h + 1.0);
                Eigen::VectorXd f_init(n*n);

                for (int i = 0; i < n*n; i++)
		        {
                    double x_i = (i % n)*h;
                    double y_i = (i / n)*h;

                    if (i % n == 0 || i % n == n-1 || i / n == 0 || i / n == n-1)
                    {
                        f_init(i) = u(x_i, y_i);
                    }
                    else
                    {
                        f_init(i) = u.Laplace(x_i, y_i);
                    }
                }

                return f_init;
            }
            if (boundary == "Neumann")
            {
                int n = static_cast<int>(1.0 / h + 1.0);
                Eigen::VectorXd f_init(n*n);

                for (int i = 0; i < n*n; i++)
		        {
                    double x_i = (i % n)*h;
                    double y_i = (i / n)*h;

                    if (i == 0 || i == n-1 || i == n*n-n || i == n*n-1)
                    {
                        f_init(i) = u(x_i, y_i)/(h*h);
                    }
                    else
                    {
                        if (i % n == 0 || i % n == n-1 || i / n == 0 || i / n == n-1)
                        {
                            if (i % n == 0) // Left Boundary
                            {
                                f_init(i) = u.diff_x(x_i, y_i)/h;
                            }
                            if (i % n == n-1) // Right Boundary
                            {
                                f_init(i) = u.diff_x(x_i, y_i)/h;
                            }
                            if (i / n == 0) // Lower Boundary
                            {
                                if (i == 1)
                                {
                                    f_init(i) = u(h, 0.0)/(h*h);
                                }
                                else
                                {
                                    f_init(i) = u.diff_y(x_i, y_i)/h;
                                } 
                            }
                            if (i / n == n-1) // Upper Boundary
                            {
                                f_init(i) = u.diff_y(x_i, y_i)/h;
                            }
                        }
                        else
                        {
                            f_init(i) = u.Laplace(x_i, y_i);
                        }
                    }
		        }

                return f_init;
            }
            if (boundary == "Mixed")
            {
                int n = static_cast<int>(1.0 / h + 1.0);
                Eigen::VectorXd f_init(n*n);

                for (int i = 0; i < n*n; i++)
		        {
                    double x_i = (i % n)*h;
                    double y_i = (i / n)*h;

                    if (i == 0 || i == n-1 || i == n*n-n || i == n*n-1)
                    {
                        f_init(i) = u(x_i, y_i);
                    }
                    else
                    {
                        if (i % n == 0 || i % n == n-1 || i / n == 0 || i / n == n-1)
                        {
                            if (i % n == 0) // Left Boundary
                            {
                                f_init(i) = u.diff_x(x_i, y_i);
                            }
                            if (i % n == n-1) // Right Boundary
                            {
                                f_init(i) = u.diff_x(x_i, y_i);
                            }
                            if (i / n == 0) // Lower Boundary
                            {
                                f_init(i) = u(x_i, y_i);
                            }
                            if (i / n == n-1) // Upper Boundary
                            {
                                f_init(i) = u(x_i, y_i);
                            }
                        }
                        else
                        {
                            f_init(i) = u.Laplace(x_i, y_i);
                        }
                    }
		        }

                return f_init;
            }
        }
        Eigen::VectorXd vec;
        return vec;
    }
};


/*
    Solver
*/
template <int dim>
class Solver
{
private:
    int v1, v2;
    Function<dim> &u;          // Grid function
    std::string boundary;      // Boundary conditions
    std::string resop;         // restriction operator
    std::string interop;       // interpolation operator
    int maxiter;               // Maximum number of iterations
    double relacc;             // Relative convergence precision
    double h;                  // Grid width
    std::vector<double> initv; // Initial vector for iteration
    Eigen::VectorXd v;         // Compute solution
    std::vector<Eigen::SparseMatrix<double>> A; 
    std::vector<double> Re;    // residual
    double CPU_time;           // Timing

    // V-cycle
    Eigen::VectorXd solve_Vcycle(Eigen::VectorXd v_h, Eigen::VectorXd f_h, double H)
    {
        Eigen::SparseMatrix<double> A_h = A[log2(static_cast<int>(H/h))];
        Weighted_Jacobi<dim> J1(A_h, f_h, v_h, v1, 2.0/3.0);
        v_h = J1.solve();
        if (H >= 0.5)
        {
            Weighted_Jacobi<dim> J2(A_h, f_h, v_h, v2, 2.0/3.0);
            return J2.solve();
        }
        else
        {
            Transformation<dim> T1(f_h - A_h * v_h, H);
            Eigen::VectorXd f_2h = T1.Restriction(resop, boundary);
            Eigen::VectorXd v_2h = Eigen::VectorXd::Zero(f_2h.size());
            v_2h = solve_Vcycle(v_2h, f_2h, 2 * H);
            Transformation<dim> T2(v_2h, 2 * H);
            v_h = v_h + T2.Interpolation(interop, boundary);
            Weighted_Jacobi<dim> J2(A_h, f_h, v_h, v2, 2.0/3.0);
            return J2.solve();
        }
    }

    Eigen::VectorXd solve_FMG(Eigen::VectorXd f_h, double H)
    {
        Eigen::VectorXd v_h;
        if (H >= 0.5)
        {
            v_h = Eigen::VectorXd::Zero(f_h.size());
            v_h = solve_Vcycle(v_h, f_h, H);
            return v_h;
        }
        else
        {
            Transformation<dim> T1(f_h, H);
            Eigen::VectorXd f_2h = T1.Restriction(resop, boundary);
            Eigen::VectorXd v_2h = Eigen::VectorXd::Zero(f_2h.size());
            v_2h = solve_FMG(f_2h, 2 * H);
            Transformation<dim> T2(v_2h, 2 * H);
            v_h = T2.Interpolation(interop, boundary);
            v_h = solve_Vcycle(v_h, f_h, H);

            return v_h;
        }
        
    }


public:
    Solver(Function<dim> &_u,
            std::string _boundary,
            std::string _resop,
            std::string _interop,
            double _maxiter,
            double _relacc,
            double _h,
            std::vector<double> _initv)
        : u(_u),
          boundary(_boundary),
          resop(_resop),
          interop(_interop),
          maxiter(_maxiter),
          relacc(_relacc),
          h(_h),
          initv(_initv) 
    {
        if (dim == 1)
        {
            v1 = 6;
            v2 = 6;
        }
        if (dim == 2)
        {
            v1 = 10;
            v2 = 10;
            if (dim == 2 && boundary == "Neumann")
            {
                v1 *= 5;
                v2 *= 5;
            }
        }

        double h_A = h;
        while (h_A <= 0.5)
        {
            Coefficient<dim> C(u, boundary, h_A);
            Eigen::SparseMatrix<double> A_h = C.get_Ah();
            A.push_back(A_h);
            h_A = h_A*2;
        }

    }

    void Vcycle()
    {
        Eigen::VectorXd init_f;
        Eigen::SparseMatrix<double> init_A;
        Coefficient<dim> C(u, boundary, h);

        init_f = C.get_initf();
        init_A = C.get_Ah();
        Eigen::VectorXd init_v = Eigen::VectorXd::Map(initv.data(), initv.size());
        v = init_v;

        int iter = 0;
        double rel_res_norm = 1.0;

        // V-cycle
        auto start_time = std::chrono::high_resolution_clock::now();
        while (iter < maxiter && rel_res_norm >= relacc)
        {
            v = solve_Vcycle(v, init_f, h);
            Eigen::VectorXd r = init_f - init_A * v;
            rel_res_norm = r.norm() / init_f.norm();
            if (dim == 1)
            {
                Re.push_back(h*r.lpNorm<1>());
            }
            if (dim == 2)
            {
                Re.push_back(h*h*r.lpNorm<1>());
            }
            iter++;
        }
        auto end_time = std::chrono::high_resolution_clock::now();
        CPU_time = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time).count();
    }

    void FMG()
    {
        Eigen::VectorXd init_f;
        Eigen::VectorXd f;
        Eigen::VectorXd v_temp;
        Eigen::SparseMatrix<double> init_A;
        Coefficient<dim> C(u, boundary, h);

        init_f = C.get_initf();
        init_A = C.get_Ah();
        v = Eigen::VectorXd::Zero(init_f.size());

        int iter = 0;
        double rel_res_norm = 1.0;
        v = solve_Vcycle(v, init_f, h);
        f = init_f - init_A*v;

        // FMG 
        auto start_time = std::chrono::high_resolution_clock::now();
        while (iter < maxiter && rel_res_norm >= relacc)
        {
            v_temp = solve_FMG(f, h);
            v = v + v_temp;
            f = f - init_A * v_temp;
            rel_res_norm = f.norm() / init_f.norm();
            if (dim == 1)
            {
                Re.push_back(sqrt(h)*f.norm());
            }
            if (dim == 2)
            {
                Re.push_back(h*f.norm());
            }
            iter++;
        }
        auto end_time = std::chrono::high_resolution_clock::now();
        CPU_time = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time).count();
    }

    std::vector<double> get_residual() // Return the one-norm of the residual of each iteration
    {
        return Re;
    }
    double get_errornorm() // Return the infinity norm of the error
    {
        Eigen::VectorXd U_hat = Eigen::VectorXd::Zero(v.size()); 

        if (dim == 1)
        {
            for (int i = 0; i < v.size(); ++i)
            {
                U_hat(i) = u(i*h);
            }
        }
        if (dim == 2)
        {
            for (int i = 0; i < v.size(); ++i)
            {
                int N = static_cast<int>(1.0 / h + 1.0); 
                double x_i = (i % N)*h;
                double y_i = (i / N)*h;
                U_hat(i) = u(x_i, y_i);
            }
        }

        double normInf = (v - U_hat).lpNorm<Eigen::Infinity>();

        return normInf;

    }
    double get_CPUtime() // Return the timing in milliseconds
    {
        return CPU_time;
    }
};

