#include <iostream>
#include <fstream>
#include "MC.h"
#include <cmath>
#include <algorithm>
#include <vector>
#include <jsoncpp/json/json.h>

#define DIM 1
#define FUNCTION Fun1

template <int dim>
class Fun1 : public Function<dim>
{
public:
    double operator()(double x) override
    {
        return exp(sin(x));
    }
    double diff(double x) override
    {
        return exp(sin(x))*cos(x);
    }
    double Laplace(double x) override
    {
        return exp(sin(x))*sin(x)-exp(sin(x))*cos(x)*cos(x);
    }

    double operator()(double _x, double _y) override
    {
        return exp(_y+sin(_x));
    }
    double diff_x(double _x, double _y) override 
    {
        return cos(_x)*exp(_y+sin(_x));  
    }
    double diff_y(double _x, double _y) override
     {
        return exp(_y+sin(_x));  
    }
    double Laplace(double _x, double _y) override
    {
        return sin(_x)*exp(_y+sin(_x)) - cos(_x)*cos(_x)*exp(_y+sin(_x)) - exp(_y+sin(_x));
    }
};

template <int dim>
class Fun2 : public Function<dim>
{
public:
    double operator()(double x) override
    {
        return sin(x)+1;
    }
    double diff(double x) override
    {
        return cos(x);
    }
    double Laplace(double x) override
    {
        return sin(x);
    }

    double operator()(double _x, double _y) override
    {
        return sin(_x + _y)+1;
    }
    double diff_x(double _x, double _y) override 
    {
        return cos(_x + _y);  
    }
    double diff_y(double _x, double _y) override
     {
        return cos(_x + _y);  
    }
    double Laplace(double _x, double _y) override
    {
        return 2.0*sin(_x + _y);
    }
};

template <int dim>
class Fun3 : public Function<dim>
{
public:
    double operator()(double x) override
    {
        return exp(-x);
    }
    double diff(double x) override
    {
        return -exp(-x);
    }
    double Laplace(double x) override
    {
        return -exp(-x);
    }

    double operator()(double _x, double _y) override
    {
        return exp(-_x+_y);
    }
    double diff_x(double _x, double _y) override
    {
        return -exp(-_x+_y);  
    }
    double diff_y(double _x, double _y) override
    {
        return exp(-_x+_y);  
    }
    double Laplace(double _x, double _y) override
    {
        return -2.0*exp(-_x+_y);
    }
};

int main()
{
    Json::Reader reader;
	Json::Value root;
 
	std::ifstream in("test.json", std::ios::binary);
	if (!in.is_open())
	{
		std::cout << "error: cannot open file." << std::endl;
		return -1;
    }

    int maximum_iteration_number;
    double relative_accuracy;
    std::vector<int> n;
    std::string boundary_condition, restriction_operators, interpolation_operators, cycles, initial_guess;

    if (reader.parse(in, root))
	{
        for (int i = 0; i < root["n"].size(); i++)
	    {
	        int temp_n = root["n"][i].asInt();
	        n.push_back(temp_n);
	    }
        boundary_condition = root["boundary_condition"].asString();
        restriction_operators = root["restriction_operator"].asString();
        interpolation_operators = root["interpolation_operator"].asString();
        cycles = root["cycle"].asString();
        initial_guess = root["initial_guess"].asString();
        maximum_iteration_number = root["maximum_iteration_number"].asInt();
        relative_accuracy = root["relative_accuracy"].asDouble();
    }

    const int dim = DIM;
    FUNCTION<dim> f;
    std::vector<double> Error_norm, Time_cost, Point1_error, Point2_error, Point3_error, Point4_error;
    for (int i = 0; i < n.size(); i++)
	{
        std::vector<double> v;
        if (dim == 1 && initial_guess == "null vector") 
        {
            v = std::vector<double>(n[i]+1, 0.0);
        } 
        if (dim == 2 && initial_guess == "null vector") 
        {
            v = std::vector<double>((n[i]+1)*(n[i]+1), 0.0);
        } 
        Solver<dim> V(f,boundary_condition, restriction_operators, interpolation_operators, maximum_iteration_number, relative_accuracy, 1.0/n[i], v);
        if (cycles == "V")
        {
            V.Vcycle();
        }
        if (cycles == "FMG")
        {
            V.FMG();
        }

        std::vector<double> R = V.get_residual();
        double error_norm = V.get_errornorm();
        //std::vector<double> point_error = V.get_4points_result();
        double time_cost = V.get_CPUtime();

        Error_norm.push_back(error_norm);
        Time_cost.push_back(time_cost);
        R.erase(R.begin());

        std::cout << "-------------------n = " << n[i] << "-------------------" << std::endl;
        std::cout << "Iteration number: " << R.size() << std::endl;
        std::cout << "Residual's norm AND relative reduction rate between each cycle: " << std::endl;
        std::cout << R[0] << std::endl;
        for (int j = 1; j < R.size(); ++j)
        {
            std::cout << R[j] << " & " << round(R[j]/R[j-1]*1000)/1000 << std::endl;
        }
    }

    std::cout << "------------------errornorm report-------------------" << std::endl;
    std::cout << "Errornorm of each n: " << std::endl;
    for (int i = 0; i < Error_norm.size(); ++i)
    {
        std::cout << Error_norm[i] << ", ";
    }
    std::cout << std::endl;
    std::cout << "Relative convergence rate(1og 2) between each n: " << std::endl;
    for (int i = 1; i < Error_norm.size(); ++i)
    {
        std::cout << log2(Error_norm[i-1]/Error_norm[i]) << ", ";
    }
    std::cout << std::endl;

    std::cout << "------------------CPU time report-------------------" << std::endl;
    std::cout << "CPU time of each n(ms): " << std::endl;
    for (int i = 0; i < Time_cost.size(); ++i)
    {
        std::cout << Time_cost[i] << ", ";
    }
    std::cout << std::endl;

    return 0;
}

