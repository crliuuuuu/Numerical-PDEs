#include <iostream>
#include <fstream>
#include "FD.h"
#include <cmath>
#include <algorithm>
#include <chrono>
#include <vector>
#include <jsoncpp/json/json.h>


class Fun1 : public Function
{
public:
    double operator()(double _x, double _y)
    {
        return std::exp(_y+std::sin(_x));
    }
    double diff_x(double _x, double _y) 
    {
        return std::cos(_x)*std::exp(_y+std::sin(_x));  
    }
    double diff_y(double _x, double _y) 
    {
        return std::exp(_y+std::sin(_x));  
    }
    double Laplace(double _x, double _y)
    {
        return std::sin(_x)*std::exp(_y+std::sin(_x)) - std::cos(_x)*std::cos(_x)*std::exp(_y+std::sin(_x)) - std::exp(_y+std::sin(_x));
    }
};

class Fun2 : public Function
{
public:
    double operator()(double _x, double _y)
    {
        return std::exp(-_x*_y);
    }
    double diff_x(double _x, double _y) 
    {
        return std::exp(-_x*_y)*(-_y);  
    }
    double diff_y(double _x, double _y) 
    {
        return std::exp(-_x*_y)*(-_x);  
    }
    double Laplace(double _x, double _y)
    {
        return -std::exp(-_x*_y)*_x*_x-std::exp(-_x*_y)*_y*_y;
    }
};

class Fun3 : public Function
{
public:
    double operator()(double _x, double _y)
    {
        return std::sin(_x*_y);
    }
    double diff_x(double _x, double _y) 
    {
        return std::cos(_x*_y)*_y;  
    }
    double diff_y(double _x, double _y) 
    {
        return std::cos(_x*_y)*_x;  
    }
    double Laplace(double _x, double _y)
    {
        return std::sin(_x*_y)*_y*_y + std::sin(_x*_y)*_x*_x;
    }
};

int main()
{
    Json::Reader reader;
	Json::Value root;
    std::ofstream fout;
 
	std::ifstream in("test.json", std::ios::binary);
	if (!in.is_open())
	{
		std::cout << "error: cannot open file." << std::endl;
		return -1;
    }

    std::vector<double> h, circle;
    int function_label, boundary_condition;

    if (reader.parse(in, root))
	{
        function_label = root["function_label"].asInt();
        boundary_condition = root["boundary_condition"].asInt();
        for (int i = 0; i < root["n"].size(); i++)
	    {
	        int temp_n = root["n"][i].asInt();
	        h.push_back(1.0/temp_n);
	    }
        for (int i = 0; i < root["circle"].size(); i++)
	    {
	        double temp_n = root["circle"][i].asDouble();
	        circle.push_back(temp_n);
	    }
    }

    fout.open("test_results.txt");
    fout << "Results of Function" << function_label << ", condition " << boundary_condition << ": " << std::endl;

    for (int i = 0; i < h.size(); i++)
	{
        if (circle[2] == 0.0)
        {
            if (function_label == 1)
            {
                Fun1 F;
                FD_regular FD(F, h[i], boundary_condition);
                auto start_time = std::chrono::high_resolution_clock::now();
                FD.solve();
                auto end_time = std::chrono::high_resolution_clock::now();
                double CPU_time = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time).count();
                std::cout << CPU_time << " ";
                std::vector<double> error = FD.get_errornorm();
                fout << "n = " << static_cast<int>(1.0 / h[i]) << std::endl;
                fout << "error norms(1, 2 and infinity) are " << error[0] << ", " << error[1] << ", " << error[2] << std::endl;
                fout << "results(from left to right, from bottom to top): " << std::endl;
                for (int j = 0; j <= static_cast<int>(1.0 / h[i]); j++)
	            {
                    for (int k = 0; k <= static_cast<int>(1.0 / h[i]); k++)
	                {
                        if (k == 0)
                        {
                            fout << "[" ;
                        }
                        fout << FD.get_result(k*h[i], j*h[i]);
                        if (k == static_cast<int>(1.0 / h[i]))
                        {
                            fout << "], " << std::endl;
                        }
                        else
                        {
                            fout << ", ";
                        }
                    } 
                }
                fout << "error at (0.125,0.125),(0.875,0.875),(0.125,0.875),(0.875,0.125): " << std::endl;
                fout << std::abs(FD.get_result(0.125,0.125) - F(0.125,0.125)) << ", ";
                fout << std::abs(FD.get_result(0.875,0.875) - F(0.875,0.875)) << ", ";
                fout << std::abs(FD.get_result(0.125,0.875) - F(0.125,0.875)) << ", ";
                fout << std::abs(FD.get_result(0.875,0.125) - F(0.875,0.125)) << std::endl;
            }
            if (function_label == 2)
            {
                Fun2 F;
                FD_regular FD(F, h[i], boundary_condition);
                FD.solve();
                std::vector<double> error = FD.get_errornorm();
                fout << "n = " << static_cast<int>(1.0 / h[i]) << std::endl;
                fout << "error norms(1, 2 and infinity) are " << error[0] << ", " << error[1] << ", " << error[2] << std::endl;
                fout << "results(from left to right, from bottom to top): " << std::endl;
                for (int j = 0; j <= static_cast<int>(1.0 / h[i]); j++)
	            {
                    for (int k = 0; k <= static_cast<int>(1.0 / h[i]); k++)
	                {
                        if (k == 0)
                        {
                            fout << "[" ;
                        }
                        fout << FD.get_result(k*h[i], j*h[i]);
                        if (k == static_cast<int>(1.0 / h[i]))
                        {
                            fout << "], " << std::endl;
                        }
                        else
                        {
                            fout << ", ";
                        }
                    } 
                }
                fout << "error at (0.125,0.125),(0.875,0.875),(0.125,0.875),(0.875,0.125): " << std::endl;
                fout << std::abs(FD.get_result(0.125,0.125) - F(0.125,0.125)) << ", ";
                fout << std::abs(FD.get_result(0.875,0.875) - F(0.875,0.875)) << ", ";
                fout << std::abs(FD.get_result(0.125,0.875) - F(0.125,0.875)) << ", ";
                fout << std::abs(FD.get_result(0.875,0.125) - F(0.875,0.125)) << std::endl;
            }
            if (function_label == 3)
            {
                Fun3 F;
                FD_regular FD(F, h[i], boundary_condition);
                FD.solve();
                std::vector<double> error = FD.get_errornorm();
                fout << "n = " << static_cast<int>(1.0 / h[i]) << std::endl;
                fout << "error norms(1, 2 and infinity) are " << error[0] << ", " << error[1] << ", " << error[2] << std::endl;
                fout << "results(from left to right, from bottom to top): " << std::endl;
                for (int j = 0; j <= static_cast<int>(1.0 / h[i]); j++)
	            {
                    for (int k = 0; k <= static_cast<int>(1.0 / h[i]); k++)
	                {
                        if (k == 0)
                        {
                            fout << "[" ;
                        }
                        fout << FD.get_result(k*h[i], j*h[i]);
                        if (k == static_cast<int>(1.0 / h[i]))
                        {
                            fout << "], " << std::endl;
                        }
                        else
                        {
                            fout << ", ";
                        }
                    } 
                }
                fout << "error at (0.125,0.125),(0.875,0.875),(0.125,0.875),(0.875,0.125): " << std::endl;
                fout << std::abs(FD.get_result(0.125,0.125) - F(0.125,0.125)) << ", ";
                fout << std::abs(FD.get_result(0.875,0.875) - F(0.875,0.875)) << ", ";
                fout << std::abs(FD.get_result(0.125,0.875) - F(0.125,0.875)) << ", ";
                fout << std::abs(FD.get_result(0.875,0.125) - F(0.875,0.125)) << std::endl;
            }
        }
        else
        {
            if (function_label == 1)
            {
                Fun1 F;
                Circle D({circle[0],circle[1]}, circle[2]);
                FD_irregular FDI(F, h[i], boundary_condition, D);
                FDI.solve();
                std::vector<double> error = FDI.get_errornorm();
                fout << "n = " << static_cast<int>(1.0 / h[i]) << std::endl;
                fout << "error norms(1, 2 and infinity) are " << error[0] << ", " << error[1] << ", " << error[2] << std::endl;
                fout << "error at (0.125,0.125),(0.875,0.875),(0.125,0.875),(0.875,0.125): " << std::endl;
                fout << std::abs(FDI.get_result(0.125,0.125) - F(0.125,0.125)) << ", ";
                fout << std::abs(FDI.get_result(0.875,0.875) - F(0.875,0.875)) << ", ";
                fout << std::abs(FDI.get_result(0.125,0.875) - F(0.125,0.875)) << ", ";
                fout << std::abs(FDI.get_result(0.875,0.125) - F(0.875,0.125)) << std::endl;
            }
            if (function_label == 2)
            {
                Fun2 F;
                Circle D({circle[0],circle[1]}, circle[2]);
                FD_irregular FDI(F, h[i], boundary_condition, D);
                FDI.solve();
                std::vector<double> error = FDI.get_errornorm();
                fout << "n = " << static_cast<int>(1.0 / h[i]) << std::endl;
                fout << "error norms(1, 2 and infinity) are " << error[0] << ", " << error[1] << ", " << error[2] << std::endl;
                fout << "error at (0.125,0.125),(0.875,0.875),(0.125,0.875),(0.875,0.125): " << std::endl;
                fout << std::abs(FDI.get_result(0.125,0.125) - F(0.125,0.125)) << ", ";
                fout << std::abs(FDI.get_result(0.875,0.875) - F(0.875,0.875)) << ", ";
                fout << std::abs(FDI.get_result(0.125,0.875) - F(0.125,0.875)) << ", ";
                fout << std::abs(FDI.get_result(0.875,0.125) - F(0.875,0.125)) << std::endl;
            }
            if (function_label == 3)
            {
                Fun3 F;
                Circle D({circle[0],circle[1]}, circle[2]);
                FD_irregular FDI(F, h[i], boundary_condition, D);
                FDI.solve();
                std::vector<double> error = FDI.get_errornorm();
                fout << "n = " << static_cast<int>(1.0 / h[i]) << std::endl;
                fout << "error norms(1, 2 and infinity) are " << error[0] << ", " << error[1] << ", " << error[2] << std::endl;
                fout << "error at (0.125,0.125),(0.875,0.875),(0.125,0.875),(0.875,0.125): " << std::endl;
                fout << std::abs(FDI.get_result(0.125,0.125) - F(0.125,0.125)) << ", ";
                fout << std::abs(FDI.get_result(0.875,0.875) - F(0.875,0.875)) << ", ";
                fout << std::abs(FDI.get_result(0.125,0.875) - F(0.125,0.875)) << ", ";
                fout << std::abs(FDI.get_result(0.875,0.125) - F(0.875,0.125)) << std::endl;
            }
        }
    }

    fout.close();

    return 0;
}

