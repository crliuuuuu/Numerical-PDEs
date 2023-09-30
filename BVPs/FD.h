#include <iostream>
#include <cmath>
#include <algorithm>
#include <vector>
#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/Sparse>
#include <eigen3/Eigen/Core>


class Function
{
public:
    virtual double operator()(double _x, double _y) = 0; 
    virtual double diff_x(double _x, double _y) // Return the value of the partial derivative of function u with respect to x at (_x,_y)
    {
        return 0;  
    }
    virtual double diff_y(double _x, double _y) // Return the value of the partial derivative of function u with respect to y at (_x,_y)
    {
        return 0;  
    }
    virtual double Laplace(double _x, double _y) // Return the value of function u at (_x,_y) after the -Laplace operator is applied
    {
        return 0;  
    }
};

/*
    Define FD method class, including the use of FD (Finite Difference) method on the unit square and the unit square excluding a circle
*/
class FD
{
public:
    virtual void solve() = 0;  
    virtual double get_result(double _x, double _y) // Return the calculated value at grid point (_x,_y)
    {
        return 0.0;  
    }  
    virtual std::vector<double> get_errornorm() // Return the 1, 2, and infinity norms of error
    {
        return std::vector<double>();  
    }  
};

class Circle 
{
private:
    std::vector<double> center_;  
    double radius_;              
public:
    Circle(std::vector<double> center, double radius) : center_(center), radius_(radius) {}

    // Calculate the unit outward vector from point P to the center of the circle
    std::vector<double> get_normal(std::vector<double> P) 
    {
        if (P.size() != 2) 
        {
            std::cerr<< "Invalid dimensions for point or circle." <<std::endl;
		    exit(-1);
        }
        std::vector<double> normal = {P[0] - center_[0], P[1] - center_[1]};
        double length = std::sqrt(std::pow(normal[0], 2) + std::pow(normal[1], 2));
        normal[0] /= length;
        normal[1] /= length;
        return normal;
    }

    // Calculate the two intersection points of the line and the circle
    std::vector<double> get_intersection(bool is_vertical, double value) 
    {
        if (is_vertical)  
        {
            if (std::abs(center_[0] - value) >= radius_) {
                return {-1.0, -1.0};
            }
            double x = value;
            double y1 = center_[1] + std::sqrt(std::pow(radius_, 2) - std::pow(x - center_[0], 2));
            double y2 = center_[1] - std::sqrt(std::pow(radius_, 2) - std::pow(x - center_[0], 2));
            return {y1, y2};
        } else 
        {
            if (std::abs(center_[1] - value) >= radius_) {
                return {-1.0, -1.0};
            }
            double y = value;
            double x1 = center_[0] + std::sqrt(std::pow(radius_, 2) - std::pow(y - center_[1], 2));
            double x2 = center_[0] + std::sqrt(std::pow(radius_, 2) - std::pow(y - center_[1], 2));
            return {x1, x2};
        }
    }

    // Return the center of the circle
    std::vector<double> get_center()
    {
        return center_ ;
    }

    // Return the radius of the circle
    double get_radius()
    {
        return radius_ ;
    }
};

/*
    Solve the differential equation -laplace u = f on (0,1)*(0,1) using the FD method, where f and boundary conditions are given
*/
class FD_regular : public FD
{
private:
    Function &u;
    double h;
    int condition, N; // N: Number of grid points on each row or column (including boundary)
    Eigen::VectorXd U; // U: The computed solution
public:
    /*
    u: The true function value, used to calculate f and boundary conditions, and used for error analysis compared with the computed value
    h: Grid width
    condition: Boundary conditions (1 for Dirichlet, 2 for Neumann, 3 for mixed)
    */
    FD_regular(Function &_u, double _h, int _condition): u(_u), h(_h), condition(_condition)
    {
        if ((condition != 1) && (condition != 2) && (condition != 3))
	    {
		    std::cerr<< "Invalid condition." <<std::endl;
		    exit(-1);
	    }

        N = static_cast<int>(1.0 / h + 1.0); 
    }

    // solve AU = F
    void solve()
    {
        if (condition == 1) // Dirichlet
        {
            Eigen::SparseMatrix<double> A(N*N, N*N);
            std::vector<Eigen::Triplet<double> > tripletlist; 
            Eigen::MatrixXd F(N*N,1); 
            F = Eigen::MatrixXd::Zero(N*N, 1); 

            for (int i = 0; i < N*N; i++)
		    {
                double x_i = (i % N)*h;
                double y_i = (i / N)*h;

                if (i % N == 0 || i % N == N-1 || i / N == 0 || i / N == N-1)
                {
                    tripletlist.push_back(Eigen::Triplet<double>(i, i, 1.0));
                    F(i, 0) = u(x_i, y_i);
                }
                else
                {
                    tripletlist.push_back(Eigen::Triplet<double>(i, i, 4.0/(h*h)));
                    tripletlist.push_back(Eigen::Triplet<double>(i, i-1, -1.0/(h*h)));
                    tripletlist.push_back(Eigen::Triplet<double>(i, i+1, -1.0/(h*h)));
                    tripletlist.push_back(Eigen::Triplet<double>(i, i-N, -1.0/(h*h)));
                    tripletlist.push_back(Eigen::Triplet<double>(i, i+N, -1.0/(h*h)));
                    F(i, 0) = u.Laplace(x_i, y_i);
                }
		    }

            A.setFromTriplets(tripletlist.begin(), tripletlist.end());
            A.makeCompressed();
            Eigen::SparseLU<Eigen::SparseMatrix<double> > Solver_sparse;
            Solver_sparse.compute(A);
            U = Solver_sparse.solve(F);
        }

        if (condition == 2) // Neumann
        {
            Eigen::SparseMatrix<double> A(N*N, N*N); 
            std::vector<Eigen::Triplet<double> > tripletlist; 
            Eigen::MatrixXd F(N*N,1); 
            F = Eigen::MatrixXd::Zero(N*N, 1); 

            for (int i = 0; i < N*N; i++)
		    {
                double x_i = (i % N)*h;
                double y_i = (i / N)*h;

                if (i == 0 || i == N-1 || i == N*N-N || i == N*N-1)
                {
                    tripletlist.push_back(Eigen::Triplet<double>(i, i, 1.0));
                    F(i, 0) = u(x_i, y_i);
                }
                else
                {
                    if (i % N == 0 || i % N == N-1 || i / N == 0 || i / N == N-1)
                    {
                        if (i % N == 0) // Left boundary
                        {
                            tripletlist.push_back(Eigen::Triplet<double>(i, i, -3.0/(2.0*h)));
                            tripletlist.push_back(Eigen::Triplet<double>(i, i+1, 2.0/h));
                            tripletlist.push_back(Eigen::Triplet<double>(i, i+2, -1.0/(2.0*h)));
                            F(i, 0) = u.diff_x(x_i, y_i);
                        }
                        if (i % N == N-1) // Right boundary
                        {
                            tripletlist.push_back(Eigen::Triplet<double>(i, i, 3.0/(2.0*h)));
                            tripletlist.push_back(Eigen::Triplet<double>(i, i-1, -2.0/h));
                            tripletlist.push_back(Eigen::Triplet<double>(i, i-2, 1.0/(2.0*h)));
                            F(i, 0) = u.diff_x(x_i, y_i);
                        }
                        if (i / N == 0) // Lower boundary
                        {
                            if (i == 1)
                            {
                                tripletlist.push_back(Eigen::Triplet<double>(i, i, 1.0));
                                F(i, 0) = u(h, 0.0);
                            }
                            else
                            {
                                tripletlist.push_back(Eigen::Triplet<double>(i, i, -3.0/(2.0*h)));
                                tripletlist.push_back(Eigen::Triplet<double>(i, i+N, 2.0/h));
                                tripletlist.push_back(Eigen::Triplet<double>(i, i+2*N, -1.0/(2.0*h)));
                                F(i, 0) = u.diff_y(x_i, y_i);
                            }
                        }
                        if (i / N == N-1) // Upper boundary
                        {
                            tripletlist.push_back(Eigen::Triplet<double>(i, i, 3.0/(2.0*h)));
                            tripletlist.push_back(Eigen::Triplet<double>(i, i-N, -2.0/h));
                            tripletlist.push_back(Eigen::Triplet<double>(i, i-2*N, 1.0/(2.0*h)));
                            F(i, 0) = u.diff_y(x_i, y_i);
                        }
                    }
                    else
                    {
                        tripletlist.push_back(Eigen::Triplet<double>(i, i, 4.0/(h*h)));
                        tripletlist.push_back(Eigen::Triplet<double>(i, i-1, -1.0/(h*h)));
                        tripletlist.push_back(Eigen::Triplet<double>(i, i+1, -1.0/(h*h)));
                        tripletlist.push_back(Eigen::Triplet<double>(i, i-N, -1.0/(h*h)));
                        tripletlist.push_back(Eigen::Triplet<double>(i, i+N, -1.0/(h*h)));
                        F(i, 0) = u.Laplace(x_i, y_i);
                    }
                }
		    }

            A.setFromTriplets(tripletlist.begin(), tripletlist.end());
            A.makeCompressed();
            Eigen::SparseLU<Eigen::SparseMatrix<double> > Solver_sparse;
            Solver_sparse.compute(A);
            U = Solver_sparse.solve(F);
        }

        if (condition == 3) // // Mixed, e.g. left and right are Neumann, top and bottom are Dirichlet
        {
            Eigen::SparseMatrix<double> A(N*N, N*N); 
            std::vector<Eigen::Triplet<double> > tripletlist; 
            Eigen::MatrixXd F(N*N,1); 
            F = Eigen::MatrixXd::Zero(N*N, 1); 

            for (int i = 0; i < N*N; i++)
		    {

                double x_i = (i % N)*h;
                double y_i = (i / N)*h;

                if (i == 0 || i == N-1 || i == N*N-N || i == N*N-1)
                {
                    tripletlist.push_back(Eigen::Triplet<double>(i, i, 1.0));
                    F(i, 0) = u(x_i, y_i);
                }
                else
                {
                    if (i % N == 0 || i % N == N-1 || i / N == 0 || i / N == N-1)
                    {
                        if (i % N == 0) // Left boundary
                        {
                            tripletlist.push_back(Eigen::Triplet<double>(i, i, -3.0/(2.0*h)));
                            tripletlist.push_back(Eigen::Triplet<double>(i, i+1, 2.0/h));
                            tripletlist.push_back(Eigen::Triplet<double>(i, i+2, -1.0/(2.0*h)));
                            F(i, 0) = u.diff_x(x_i, y_i);
                        }
                        if (i % N == N-1) // Right boundary
                        {
                            tripletlist.push_back(Eigen::Triplet<double>(i, i, 3.0/(2.0*h)));
                            tripletlist.push_back(Eigen::Triplet<double>(i, i-1, -2.0/h));
                            tripletlist.push_back(Eigen::Triplet<double>(i, i-2, 1.0/(2.0*h)));
                            F(i, 0) = u.diff_x(x_i, y_i);
                        }
                        if (i / N == 0) // Lower boundary
                        {
                            tripletlist.push_back(Eigen::Triplet<double>(i, i, 1.0));
                            F(i, 0) = u(x_i, y_i);
                        }
                        if (i / N == N-1) // Upper boundary
                        {
                            tripletlist.push_back(Eigen::Triplet<double>(i, i, 1.0));
                            F(i, 0) = u(x_i, y_i);
                        }
                    }
                    else
                    {
                        tripletlist.push_back(Eigen::Triplet<double>(i, i, 4.0/(h*h)));
                        tripletlist.push_back(Eigen::Triplet<double>(i, i-1, -1.0/(h*h)));
                        tripletlist.push_back(Eigen::Triplet<double>(i, i+1, -1.0/(h*h)));
                        tripletlist.push_back(Eigen::Triplet<double>(i, i-N, -1.0/(h*h)));
                        tripletlist.push_back(Eigen::Triplet<double>(i, i+N, -1.0/(h*h)));
                        F(i, 0) = u.Laplace(x_i, y_i);
                    }
                }
		    }

            A.setFromTriplets(tripletlist.begin(), tripletlist.end());
            A.makeCompressed();
            Eigen::SparseLU<Eigen::SparseMatrix<double> > Solver_sparse;
            Solver_sparse.compute(A);
            U = Solver_sparse.solve(F);
        }
    }

    double get_result(double _x, double _y) 
    {
        if (std::abs(_x/h - std::round(_x/h)) > 1e-9 || std::abs(_y/h - std::round(_y/h)) > 1e-9)
        {
            std::cerr<< "Input is NOT a grid point." <<std::endl;
		    exit(-1);
        }

        return U[std::round(_y/h)*N+std::round(_x/h)];
    }  

    std::vector<double> get_errornorm() // Return the 1, 2, and infinity norms of error
    {
        Eigen::VectorXd U_hat(N*N); 
        Eigen::VectorXd E(N*N);
  
        for (int i = 0; i < N*N; ++i)
        {
            double x_i = (i % N)*h;
            double y_i = (i / N)*h;
            U_hat(i) = u(x_i, y_i);
        }
        double norm1 = h*h*(U - U_hat).lpNorm<1>();
        double norm2 = h*(U - U_hat).norm();
        double normInf = (U - U_hat).lpNorm<Eigen::Infinity>();

        return {norm1, norm2, normInf};
    }  
};

/*
    Solve the differential equation -laplace u = f on (0,1)*(0,1) excluding a circle using the FD method, where f and boundary conditions are given
*/
class FD_irregular : public FD
{
private:
    Circle c;
    Function &u;
    double h, r, cx, cy;
    int condition, N, Neq; // N: Number of grid points on each row or column (including boundary); Neq: Number of equation systems to be actually solved
    Eigen::VectorXd U; // U: The computed solution
    std::vector<std::vector<double>> L; // L: All the points to be solved and their types (see Label function)
public:
    /*
    u: The true function value, used to calculate f and boundary conditions, and used for error analysis compared with the computed value
    h: Grid width
    condition: Boundary conditions (1 for Dirichlet, 2 for Neumann, 3 for mixed)
    c: A circle cut out from the region
    */ 
    FD_irregular(Function &_u, double _h, int _condition, Circle _c): u(_u), h(_h), condition(_condition), c(_c)
    {
        if ((condition != 1) && (condition != 2) && (condition != 3))
	    {
		    std::cerr<< "Invalid condition." <<std::endl;
		    exit(-1);
	    }
        r = c.get_radius(); 
        cx = c.get_center()[0]; 
        cy = c.get_center()[1]; 
        if (cx - r < 0 || cx + r > 1 || cy - r < 0 || cy + r > 1) 
        {
            std::cerr<< "Circle is not completely inside the grid" <<std::endl;
		    exit(-1);
        }

        N = static_cast<int>(1.0 / h + 1.0); 
        L = Label();
    }

    /*
    Return a vector, returning all the points to be solved and their types.
    The i-th element is {x_i, y_i, type_i},
    where (x_i, y_i) are the coordinates of the point with sequence number i, and type_i is the type of the point:
    type_i=1 for grid points inside the region, type_i=2 for ghost points outside the region, type_i=3 for points on the square boundary, type_i=4 for grid points that fall exactly on the circle.
    */
    std::vector<std::vector<double>> Label()
    {
        for (int i = 0; i < N; i++) 
		{
            for (int j = 0; j < N; j++) 
		    {
                double x_ij = j*h;
                double y_ij = i*h;
                double dist = sqrt((x_ij - cx)*(x_ij - cx) + (y_ij - cy)*(y_ij - cy));
                if (i == 0 || i == N-1 || j == 0 || j == N-1)
                {
                    L.push_back({x_ij, y_ij, 3.0});
                }
                else
                {
                    if (dist == r)
                    {
                        L.push_back({x_ij, y_ij, 4.0});
                    }
                    if (dist > r)
                    {
                        L.push_back({x_ij, y_ij, 1.0});
                    }
                }
            }
        }
        for (int i = 0; i < N; i++) 
		{
            for (int j = 0; j < N; j++) 
		    {
                auto Find1 = std::find(L.begin(), L.end(), std::vector<double>{i*h, j*h, 1.0});
                auto Find3 = std::find(L.begin(), L.end(), std::vector<double>{i*h, j*h, 3.0});
                auto Find4 = std::find(L.begin(), L.end(), std::vector<double>{i*h, j*h, 4.0});
                auto find1 = std::find(L.begin(), L.end(), std::vector<double>{i*h, (j-1)*h, 1.0});
                auto find2 = std::find(L.begin(), L.end(), std::vector<double>{i*h, (j+1)*h, 1.0});
                auto find3 = std::find(L.begin(), L.end(), std::vector<double>{(i-1)*h, j*h, 1.0});
                auto find4 = std::find(L.begin(), L.end(), std::vector<double>{(i+1)*h, j*h, 1.0});
                if ((find1 != L.end() || find2 != L.end() || find3 != L.end() || find4 != L.end()) && Find1 == L.end() && Find3 == L.end() && Find4 == L.end())
                {
                    L.push_back({i*h, j*h, 2.0});
                }
            }
        }

        return L;
    }

    // Return the sequence number and Type of (_x,_y) corresponding to L
    std::vector<int> Find_Label(double _x, double _y)
    {
        auto find1 = std::find(L.begin(), L.end(), std::vector<double>{_x, _y, 1.0});
        if (find1 != L.end())
        {
            return {static_cast<int>(std::distance(L.begin(), find1)), 1};
        }
        auto find2 = std::find(L.begin(), L.end(), std::vector<double>{_x, _y, 2.0});
        if (find2 != L.end())
        {
            return {static_cast<int>(std::distance(L.begin(), find2)), 2};
        }
        auto find3 = std::find(L.begin(), L.end(), std::vector<double>{_x, _y, 3.0});
        if (find3 != L.end())
        {
            return {static_cast<int>(std::distance(L.begin(), find3)), 3};
        }
        auto find4 = std::find(L.begin(), L.end(), std::vector<double>{_x, _y, 4.0});
        if (find4 != L.end())
        {
            return {static_cast<int>(std::distance(L.begin(), find4)), 4};
        }
        return {-1, -1};
    }

    // solve AU = F
    void solve()
    {
        if (condition == 1) // Dirichlet
        {
            Neq = L.size();
            Eigen::SparseMatrix<double> A(Neq, Neq);
            std::vector<Eigen::Triplet<double> > tripletlist; 
            Eigen::MatrixXd F(Neq,1); 
            F = Eigen::MatrixXd::Zero(Neq, 1); 

            for (int i = 0; i < Neq; i++)
		    {
                double x_i = L[i][0];
                double y_i = L[i][1];
                int type = static_cast<int>(L[i][2]);

                if (type == 1)
                {
                    int Label_up = Find_Label(x_i, y_i+h)[0];
                    int Label_down = Find_Label(x_i, y_i-h)[0];
                    int Label_right = Find_Label(x_i+h, y_i)[0];
                    int Label_left = Find_Label(x_i-h, y_i)[0];

                    tripletlist.push_back(Eigen::Triplet<double>(i, i, 4.0/(h*h)));
                    tripletlist.push_back(Eigen::Triplet<double>(i, Label_up, -1.0/(h*h)));
                    tripletlist.push_back(Eigen::Triplet<double>(i, Label_down, -1.0/(h*h)));
                    tripletlist.push_back(Eigen::Triplet<double>(i, Label_right, -1.0/(h*h)));
                    tripletlist.push_back(Eigen::Triplet<double>(i, Label_left, -1.0/(h*h)));
                    F(i, 0) = u.Laplace(x_i, y_i);
                }
                if (type == 2)
                {
                    int Label_up = Find_Label(x_i, y_i+h)[0];
                    int Label_down = Find_Label(x_i, y_i-h)[0];
                    int Label_right = Find_Label(x_i+h, y_i)[0];
                    int Label_left = Find_Label(x_i-h, y_i)[0];

                    int type_up = Find_Label(x_i, y_i+h)[1];
                    int type_down = Find_Label(x_i, y_i-h)[1];
                    int type_right = Find_Label(x_i+h, y_i)[1];
                    int type_left = Find_Label(x_i-h, y_i)[1];

                    if (type_up == 1)
                    {
                        std::vector<double> I = c.get_intersection(true, x_i);
                        double y_intersect = 0.0;

                        if (I[0] >= y_i && I[0] <= y_i + h)
                        {
                            y_intersect = I[0];
                        }
                        else
                        {
                            y_intersect = I[1];
                        }

                        double u_intersect = u(x_i, y_intersect);
                        tripletlist.push_back(Eigen::Triplet<double>(i, i, 1.0/(y_intersect-y_i)));
                        tripletlist.push_back(Eigen::Triplet<double>(i, Label_up, 1.0/(y_i+h-y_intersect)));
                        F(i, 0) = (1.0/(y_intersect-y_i) + 1.0/(y_i+h-y_intersect))*u_intersect;
                    }
                    else if (type_down == 1)
                    {
                        std::vector<double> I = c.get_intersection(true, x_i);
                        double y_intersect = 0.0;

                        if (I[0] <= y_i && I[0] >= y_i - h)
                        {
                            y_intersect = I[0];
                        }
                        else
                        {
                            y_intersect = I[1];
                        }

                        double u_intersect = u(x_i, y_intersect);
                        tripletlist.push_back(Eigen::Triplet<double>(i, i, 1.0/(y_i-y_intersect)));
                        tripletlist.push_back(Eigen::Triplet<double>(i, Label_down, 1.0/(y_intersect-y_i+h)));
                        F(i, 0) = (1.0/(y_i-y_intersect) + 1.0/(y_intersect-y_i+h))*u_intersect;
                    }
                    else if (type_right == 1)
                    {
                        std::vector<double> I = c.get_intersection(false, y_i);
                        double x_intersect = 0.0;

                        if (I[0] >= x_i && I[0] <= x_i + h)
                        {
                            x_intersect = I[0];
                        }
                        else
                        {
                            x_intersect = I[1];
                        }

                        double u_intersect = u(x_intersect, y_i);
                        tripletlist.push_back(Eigen::Triplet<double>(i, i, 1.0/(x_intersect-x_i)));
                        tripletlist.push_back(Eigen::Triplet<double>(i, Label_right, 1.0/(x_i+h-x_intersect)));
                        F(i, 0) = (1.0/(x_intersect-x_i) + 1.0/(x_i+h-x_intersect))*u_intersect;
                    }
                    else if (type_left == 1)
                    {
                        std::vector<double> I = c.get_intersection(false, y_i);
                        double x_intersect = 0.0;

                        if (I[0] <= x_i && I[0] >= x_i - h)
                        {
                            x_intersect = I[0];
                        }
                        else
                        {
                            x_intersect = I[1];
                        }

                        double u_intersect = u(x_intersect, y_i);
                        tripletlist.push_back(Eigen::Triplet<double>(i, i, 1.0/(x_i-x_intersect)));
                        tripletlist.push_back(Eigen::Triplet<double>(i, Label_left, 1.0/(x_intersect-x_i+h)));
                        F(i, 0) = (1.0/(x_i-x_intersect) + 1.0/(x_intersect-x_i+h))*u_intersect;
                    }
                }
                if (type == 3 || type == 4)
                {
                    tripletlist.push_back(Eigen::Triplet<double>(i, i, 1.0));
                    F(i, 0) = u(x_i, y_i);
                }
		    }

            A.setFromTriplets(tripletlist.begin(), tripletlist.end());
            A.makeCompressed();
            Eigen::SparseLU<Eigen::SparseMatrix<double> > Solver_sparse;
            Solver_sparse.compute(A);
            U = Solver_sparse.solve(F);
        }
        
        if (condition == 2) // Neumann
        {
            Neq = L.size();
            Eigen::SparseMatrix<double> A(Neq, Neq);
            std::vector<Eigen::Triplet<double> > tripletlist; 
            Eigen::MatrixXd F(Neq,1); 
            F = Eigen::MatrixXd::Zero(Neq, 1); 

            for (int i = 0; i < Neq; i++)
		    {
                double x_i = L[i][0];
                double y_i = L[i][1];
                int type = static_cast<int>(L[i][2]);

                if (i == 1) 
                {
                    tripletlist.push_back(Eigen::Triplet<double>(i, i, 1.0));
                    F(i, 0) = u(x_i, y_i);
                }
                else if (type == 1)
                {
                    int Label_up = Find_Label(x_i, y_i+h)[0];
                    int Label_down = Find_Label(x_i, y_i-h)[0];
                    int Label_right = Find_Label(x_i+h, y_i)[0];
                    int Label_left = Find_Label(x_i-h, y_i)[0];

                    tripletlist.push_back(Eigen::Triplet<double>(i, i, 4.0/(h*h)));
                    tripletlist.push_back(Eigen::Triplet<double>(i, Label_up, -1.0/(h*h)));
                    tripletlist.push_back(Eigen::Triplet<double>(i, Label_down, -1.0/(h*h)));
                    tripletlist.push_back(Eigen::Triplet<double>(i, Label_right, -1.0/(h*h)));
                    tripletlist.push_back(Eigen::Triplet<double>(i, Label_left, -1.0/(h*h)));
                    F(i, 0) = u.Laplace(x_i, y_i);
                }
                else if (type == 2 || type == 4)
                {
                    std::vector<double> direction = c.get_normal({x_i, y_i});
                    if (direction[0] <= 0 && direction[1] < 0)
                    {
                        int Label_left1 = Find_Label(x_i-h, y_i)[0];
                        int Label_left2 = Find_Label(x_i-2.0*h, y_i)[0];
                        int Label_down1 = Find_Label(x_i, y_i-h)[0];
                        int Label_down2 = Find_Label(x_i, y_i-2.0*h)[0];

                        tripletlist.push_back(Eigen::Triplet<double>(i, i, (3.0*direction[0])/(2.0*h)+(3.0*direction[1])/(2.0*h)));
                        tripletlist.push_back(Eigen::Triplet<double>(i, Label_left1, (-2.0*direction[0])/(1.0*h)));
                        tripletlist.push_back(Eigen::Triplet<double>(i, Label_left2, (1.0*direction[0])/(2.0*h)));
                        tripletlist.push_back(Eigen::Triplet<double>(i, Label_down1, (-2.0*direction[1])/(1.0*h)));
                        tripletlist.push_back(Eigen::Triplet<double>(i, Label_down2, (1.0*direction[1])/(2.0*h)));
                        F(i, 0) = u.diff_x(x_i, y_i)*direction[0] + u.diff_y(x_i, y_i)*direction[1];
                    }
                    if (direction[0] <= 0 && direction[1] >= 0)
                    {
                        int Label_left1 = Find_Label(x_i-h, y_i)[0];
                        int Label_left2 = Find_Label(x_i-2.0*h, y_i)[0];
                        int Label_up1 = Find_Label(x_i, y_i+h)[0];
                        int Label_up2 = Find_Label(x_i, y_i+2.0*h)[0];

                        tripletlist.push_back(Eigen::Triplet<double>(i, i, (3.0*direction[0])/(2.0*h)+(-3.0*direction[1])/(2.0*h)));
                        tripletlist.push_back(Eigen::Triplet<double>(i, Label_left1, (-2.0*direction[0])/(1.0*h)));
                        tripletlist.push_back(Eigen::Triplet<double>(i, Label_left2, (1.0*direction[0])/(2.0*h)));
                        tripletlist.push_back(Eigen::Triplet<double>(i, Label_up1, (2.0*direction[1])/(1.0*h)));
                        tripletlist.push_back(Eigen::Triplet<double>(i, Label_up2, (-1.0*direction[1])/(2.0*h)));
                        F(i, 0) = u.diff_x(x_i, y_i)*direction[0] + u.diff_y(x_i, y_i)*direction[1];
                    }
                    if (direction[0] > 0 && direction[1] < 0)
                    {
                        int Label_right1 = Find_Label(x_i+h, y_i)[0];
                        int Label_right2 = Find_Label(x_i+2.0*h, y_i)[0];
                        int Label_down1 = Find_Label(x_i, y_i-h)[0];
                        int Label_down2 = Find_Label(x_i, y_i-2.0*h)[0];

                        tripletlist.push_back(Eigen::Triplet<double>(i, i, (-3.0*direction[0])/(2.0*h)+(3.0*direction[1])/(2.0*h)));
                        tripletlist.push_back(Eigen::Triplet<double>(i, Label_right1, (2.0*direction[0])/(1.0*h)));
                        tripletlist.push_back(Eigen::Triplet<double>(i, Label_right2, (-1.0*direction[0])/(2.0*h)));
                        tripletlist.push_back(Eigen::Triplet<double>(i, Label_down1, (-2.0*direction[1])/(1.0*h)));
                        tripletlist.push_back(Eigen::Triplet<double>(i, Label_down2, (1.0*direction[1])/(2.0*h)));
                        F(i, 0) = u.diff_x(x_i, y_i)*direction[0] + u.diff_y(x_i, y_i)*direction[1];
                    }
                    if (direction[0] > 0 && direction[1] >= 0)
                    {
                        int Label_right1 = Find_Label(x_i+h, y_i)[0];
                        int Label_right2 = Find_Label(x_i+2.0*h, y_i)[0];
                        int Label_up1 = Find_Label(x_i, y_i+h)[0];
                        int Label_up2 = Find_Label(x_i, y_i+2.0*h)[0];

                        tripletlist.push_back(Eigen::Triplet<double>(i, i, (-3.0*direction[0])/(2.0*h)+(-3.0*direction[1])/(2.0*h)));
                        tripletlist.push_back(Eigen::Triplet<double>(i, Label_right1, (2.0*direction[0])/(1.0*h)));
                        tripletlist.push_back(Eigen::Triplet<double>(i, Label_right2, (-1.0*direction[0])/(2.0*h)));
                        tripletlist.push_back(Eigen::Triplet<double>(i, Label_up1, (2.0*direction[1])/(1.0*h)));
                        tripletlist.push_back(Eigen::Triplet<double>(i, Label_up2, (-1.0*direction[1])/(2.0*h)));
                        F(i, 0) = u.diff_x(x_i, y_i)*direction[0] + u.diff_y(x_i, y_i)*direction[1];
                    }
                }
                else if (type == 3)
                {
                    if ((std::abs(x_i) < 1e-9 || std::abs(x_i-1.0) < 1e-9) && (std::abs(y_i) < 1e-9 || std::abs(y_i-1.0) < 1e-9)) // 四个角直接赋值
                    {
                        tripletlist.push_back(Eigen::Triplet<double>(i, i, 1.0));
                        F(i, 0) = u(x_i, y_i);
                    }
                    else
                    {
                        if (std::abs(x_i) < 1e-9)
                        {
                            int Label_right1 = Find_Label(x_i+h, y_i)[0];
                            int Label_right2 = Find_Label(x_i+2.0*h, y_i)[0];
                            tripletlist.push_back(Eigen::Triplet<double>(i, i, (-3.0)/(2.0*h)));
                            tripletlist.push_back(Eigen::Triplet<double>(i, Label_right1, (2.0)/(1.0*h)));
                            tripletlist.push_back(Eigen::Triplet<double>(i, Label_right2, (-1.0)/(2.0*h)));
                            F(i, 0) = u.diff_x(x_i, y_i);
                        }
                        if (std::abs(x_i-1.0) < 1e-9)
                        {
                            int Label_left1 = Find_Label(x_i-h, y_i)[0];
                            int Label_left2 = Find_Label(x_i-2.0*h, y_i)[0];
                            tripletlist.push_back(Eigen::Triplet<double>(i, i, (3.0)/(2.0*h)));
                            tripletlist.push_back(Eigen::Triplet<double>(i, Label_left1, (-2.0)/(1.0*h)));
                            tripletlist.push_back(Eigen::Triplet<double>(i, Label_left2, (1.0)/(2.0*h)));
                            F(i, 0) = u.diff_x(x_i, y_i);
                        }
                        if (std::abs(y_i) < 1e-9)
                        {
                            int Label_up1 = Find_Label(x_i, y_i+h)[0];
                            int Label_up2 = Find_Label(x_i, y_i+2.0*h)[0];
                            tripletlist.push_back(Eigen::Triplet<double>(i, i, (-3.0)/(2.0*h)));
                            tripletlist.push_back(Eigen::Triplet<double>(i, Label_up1, (2.0)/(1.0*h)));
                            tripletlist.push_back(Eigen::Triplet<double>(i, Label_up2, (-1.0)/(2.0*h)));
                            F(i, 0) = u.diff_y(x_i, y_i);
                        }
                        if (std::abs(y_i-1.0) < 1e-9)
                        {
                            int Label_down1 = Find_Label(x_i, y_i-h)[0];
                            int Label_down2 = Find_Label(x_i, y_i-2.0*h)[0];
                            tripletlist.push_back(Eigen::Triplet<double>(i, i, (3.0)/(2.0*h)));
                            tripletlist.push_back(Eigen::Triplet<double>(i, Label_down1, (-2.0)/(1.0*h)));
                            tripletlist.push_back(Eigen::Triplet<double>(i, Label_down2, (1.0)/(2.0*h)));
                            F(i, 0) = u.diff_y(x_i, y_i);
                        }
                    }
                }
		    }

            A.setFromTriplets(tripletlist.begin(), tripletlist.end());
            A.makeCompressed();
            Eigen::SparseLU<Eigen::SparseMatrix<double> > Solver_sparse;
            Solver_sparse.compute(A);
            U = Solver_sparse.solve(F);
        }
        
        if (condition == 3) // Mixed, let's say the square boundary is Neumann, the circle boundary is Dirichlet
        {
            Neq = L.size();
            Eigen::SparseMatrix<double> A(Neq, Neq);
            std::vector<Eigen::Triplet<double> > tripletlist; 
            Eigen::MatrixXd F(Neq,1); 
            F = Eigen::MatrixXd::Zero(Neq, 1); 

            for (int i = 0; i < Neq; i++)
		    {
                double x_i = L[i][0];
                double y_i = L[i][1];
                int type = static_cast<int>(L[i][2]);

                if (type == 1)
                {
                    int Label_up = Find_Label(x_i, y_i+h)[0];
                    int Label_down = Find_Label(x_i, y_i-h)[0];
                    int Label_right = Find_Label(x_i+h, y_i)[0];
                    int Label_left = Find_Label(x_i-h, y_i)[0];

                    tripletlist.push_back(Eigen::Triplet<double>(i, i, 4.0/(h*h)));
                    tripletlist.push_back(Eigen::Triplet<double>(i, Label_up, -1.0/(h*h)));
                    tripletlist.push_back(Eigen::Triplet<double>(i, Label_down, -1.0/(h*h)));
                    tripletlist.push_back(Eigen::Triplet<double>(i, Label_right, -1.0/(h*h)));
                    tripletlist.push_back(Eigen::Triplet<double>(i, Label_left, -1.0/(h*h)));
                    F(i, 0) = u.Laplace(x_i, y_i);
                }
                else if (type == 2)
                {
                    int Label_up = Find_Label(x_i, y_i+h)[0];
                    int Label_down = Find_Label(x_i, y_i-h)[0];
                    int Label_right = Find_Label(x_i+h, y_i)[0];
                    int Label_left = Find_Label(x_i-h, y_i)[0];

                    int type_up = Find_Label(x_i, y_i+h)[1];
                    int type_down = Find_Label(x_i, y_i-h)[1];
                    int type_right = Find_Label(x_i+h, y_i)[1];
                    int type_left = Find_Label(x_i-h, y_i)[1];

                    if (type_up == 1)
                    {
                        std::vector<double> I = c.get_intersection(true, x_i);
                        double y_intersect = 0.0;

                        if (I[0] >= y_i && I[0] <= y_i + h)
                        {
                            y_intersect = I[0];
                        }
                        else
                        {
                            y_intersect = I[1];
                        }

                        double u_intersect = u(x_i, y_intersect);
                        tripletlist.push_back(Eigen::Triplet<double>(i, i, 1.0/(y_intersect-y_i)));
                        tripletlist.push_back(Eigen::Triplet<double>(i, Label_up, 1.0/(y_i+h-y_intersect)));
                        F(i, 0) = (1.0/(y_intersect-y_i) + 1.0/(y_i+h-y_intersect))*u_intersect;
                    }
                    else if (type_down == 1)
                    {
                        std::vector<double> I = c.get_intersection(true, x_i);
                        double y_intersect = 0.0;

                        if (I[0] <= y_i && I[0] >= y_i - h)
                        {
                            y_intersect = I[0];
                        }
                        else
                        {
                            y_intersect = I[1];
                        }

                        double u_intersect = u(x_i, y_intersect);
                        tripletlist.push_back(Eigen::Triplet<double>(i, i, 1.0/(y_i-y_intersect)));
                        tripletlist.push_back(Eigen::Triplet<double>(i, Label_down, 1.0/(y_intersect-y_i+h)));
                        F(i, 0) = (1.0/(y_i-y_intersect) + 1.0/(y_intersect-y_i+h))*u_intersect;
                    }
                    else if (type_right == 1)
                    {
                        std::vector<double> I = c.get_intersection(false, y_i);
                        double x_intersect = 0.0;

                        if (I[0] >= x_i && I[0] <= x_i + h)
                        {
                            x_intersect = I[0];
                        }
                        else
                        {
                            x_intersect = I[1];
                        }

                        double u_intersect = u(x_intersect, y_i);
                        tripletlist.push_back(Eigen::Triplet<double>(i, i, 1.0/(x_intersect-x_i)));
                        tripletlist.push_back(Eigen::Triplet<double>(i, Label_right, 1.0/(x_i+h-x_intersect)));
                        F(i, 0) = (1.0/(x_intersect-x_i) + 1.0/(x_i+h-x_intersect))*u_intersect;
                    }
                    else if (type_left == 1)
                    {
                        std::vector<double> I = c.get_intersection(false, y_i);
                        double x_intersect = 0.0;

                        if (I[0] <= x_i && I[0] >= x_i - h)
                        {
                            x_intersect = I[0];
                        }
                        else
                        {
                            x_intersect = I[1];
                        }

                        double u_intersect = u(x_intersect, y_i);
                        tripletlist.push_back(Eigen::Triplet<double>(i, i, 1.0/(x_i-x_intersect)));
                        tripletlist.push_back(Eigen::Triplet<double>(i, Label_left, 1.0/(x_intersect-x_i+h)));
                        F(i, 0) = (1.0/(x_i-x_intersect) + 1.0/(x_intersect-x_i+h))*u_intersect;
                    }
                }
                else if (type == 3)
                {
                    if ((std::abs(x_i) < 1e-9 || std::abs(x_i-1.0) < 1e-9) && (std::abs(y_i) < 1e-9 || std::abs(y_i-1.0) < 1e-9)) // 四个角直接赋值
                    {
                        tripletlist.push_back(Eigen::Triplet<double>(i, i, 1.0));
                        F(i, 0) = u(x_i, y_i);
                    }
                    else
                    {
                        if (std::abs(x_i) < 1e-9)
                        {
                            int Label_right1 = Find_Label(x_i+h, y_i)[0];
                            int Label_right2 = Find_Label(x_i+2.0*h, y_i)[0];
                            tripletlist.push_back(Eigen::Triplet<double>(i, i, (-3.0)/(2.0*h)));
                            tripletlist.push_back(Eigen::Triplet<double>(i, Label_right1, (2.0)/(1.0*h)));
                            tripletlist.push_back(Eigen::Triplet<double>(i, Label_right2, (-1.0)/(2.0*h)));
                            F(i, 0) = u.diff_x(x_i, y_i);
                        }
                        if (std::abs(x_i-1.0) < 1e-9)
                        {
                            int Label_left1 = Find_Label(x_i-h, y_i)[0];
                            int Label_left2 = Find_Label(x_i-2.0*h, y_i)[0];
                            tripletlist.push_back(Eigen::Triplet<double>(i, i, (3.0)/(2.0*h)));
                            tripletlist.push_back(Eigen::Triplet<double>(i, Label_left1, (-2.0)/(1.0*h)));
                            tripletlist.push_back(Eigen::Triplet<double>(i, Label_left2, (1.0)/(2.0*h)));
                            F(i, 0) = u.diff_x(x_i, y_i);
                        }
                        if (std::abs(y_i) < 1e-9)
                        {
                            int Label_up1 = Find_Label(x_i, y_i+h)[0];
                            int Label_up2 = Find_Label(x_i, y_i+2.0*h)[0];
                            tripletlist.push_back(Eigen::Triplet<double>(i, i, (-3.0)/(2.0*h)));
                            tripletlist.push_back(Eigen::Triplet<double>(i, Label_up1, (2.0)/(1.0*h)));
                            tripletlist.push_back(Eigen::Triplet<double>(i, Label_up2, (-1.0)/(2.0*h)));
                            F(i, 0) = u.diff_y(x_i, y_i);
                        }
                        if (std::abs(y_i-1.0) < 1e-9)
                        {
                            int Label_down1 = Find_Label(x_i, y_i-h)[0];
                            int Label_down2 = Find_Label(x_i, y_i-2.0*h)[0];
                            tripletlist.push_back(Eigen::Triplet<double>(i, i, (3.0)/(2.0*h)));
                            tripletlist.push_back(Eigen::Triplet<double>(i, Label_down1, (-2.0)/(1.0*h)));
                            tripletlist.push_back(Eigen::Triplet<double>(i, Label_down2, (1.0)/(2.0*h)));
                            F(i, 0) = u.diff_y(x_i, y_i);
                        }
                    }
                }
                else if (type == 4)
                {
                    tripletlist.push_back(Eigen::Triplet<double>(i, i, 1.0));
                    F(i, 0) = u(x_i, y_i);
                }
		    }

            A.setFromTriplets(tripletlist.begin(), tripletlist.end());
            A.makeCompressed();
            Eigen::SparseLU<Eigen::SparseMatrix<double> > Solver_sparse;
            Solver_sparse.compute(A);
            U = Solver_sparse.solve(F);
        }
    }

    double get_result(double _x, double _y) 
    {
        if (std::abs(_x/h - std::round(_x/h)) > 1e-9 || std::abs(_y/h - std::round(_y/h)) > 1e-9)
        {
            std::cerr<< "Input is NOT a grid point." <<std::endl;
		    exit(-1);
        }

        return U[Find_Label(_x, _y)[0]];
    }  

    std::vector<double> get_errornorm() // Return the 1, 2, and infinity norms of error
    {
        Neq = L.size();
        Eigen::VectorXd U_hat(Neq); // U_hat stores the true values
        Eigen::VectorXd E(Neq);
  
        for (int i = 0; i < Neq; ++i)
        {
            double x_i = L[i][0];
            double y_i = L[i][1];
            U_hat(i) = u(x_i, y_i);
        }
        double norm1 = h*h*(U - U_hat).lpNorm<1>();
        double norm2 = h*(U - U_hat).norm();
        double normInf = (U - U_hat).lpNorm<Eigen::Infinity>();

        return {norm1, norm2, normInf};
    }  
};

