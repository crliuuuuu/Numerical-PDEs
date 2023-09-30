#include <iostream>
#include <memory>
#include <string>
#include "TI.h"
#include <map>
#include <vector>
#include <chrono>
#include <cmath>
#include <fstream>
#include <jsoncpp/json/json.h>
#include <iomanip>

class Restricted3body_Function : public Vector_Function {
public:
    std::vector<double> operator()(std::vector<double>& input, double time) override {
        if (input.size() != 6) {
            throw std::runtime_error("Invalid input size. Expected 6 elements.");
        }

        double u1 = input[0];
        double u2 = input[1];
        double u3 = input[2];
        double u4 = input[3];
        double u5 = input[4];
        double u6 = input[5];

        double mu = 0.012277471; 
        double one_minus_mu = 1.0 - mu;

        double r1 = std::sqrt(std::pow(u2, 2) + std::pow(u3, 2) + std::pow(u1 + mu - 1, 2));
        double r2 = std::sqrt(std::pow(u2, 2) + std::pow(u3, 2) + std::pow(u1 + mu, 2));

        std::vector<double> output(6);
        output[0] = u4;
        output[1] = u5;
        output[2] = u6;
        output[3] = 2 * u5 + u1 - (mu * (u1 + mu - 1) / std::pow(r1, 3)) - (one_minus_mu * (u1 + mu) / std::pow(r2, 3));
        output[4] = -2 * u4 + u2 - (mu * u2 / std::pow(r1, 3)) - (one_minus_mu * u2 / std::pow(r2, 3));
        output[5] = -(mu * u3 / std::pow(r1, 3)) - (one_minus_mu * u3 / std::pow(r2, 3));

        return output;
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

    std::vector<double> initial_value, num_iterations;
    double period_time;
    int order;
    std::string Richardson, Method;
    bool output = false;

    if (reader.parse(in, root))
	{
        for (int i = 0; i < root["initial_value"].size(); i++)
	    {
	        double temp = root["initial_value"][i].asDouble();
	        initial_value.push_back(temp);
	    }
        for (int i = 0; i < root["one_period_step_number"].size(); i++)
	    {
	        int temp = root["one_period_step_number"][i].asInt();
	        num_iterations.push_back(temp);
	    }
        period_time = root["one_period_time"].asDouble();
        Method = root["method"].asString();
        order = root["order"].asInt();
        Richardson = root["if_Richardson"].asString();
    }

    Restricted3body_Function restricted3body_func;

    if (Richardson == "yes")
    {
        std::vector<std::vector<double>> error_h_list;
        std::vector<std::vector<double>> error_halfh_list;
        std::vector<std::vector<double>> convergence_rate_list;
        std::vector<double> average_rate;
        std::vector<double> time;
        for (int i = 0; i < num_iterations.size(); i++)
	    {
            TimeIntegratorFactory &factory = TimeIntegratorFactory::getInstance();
            double step_size = period_time / num_iterations[i];
            std::unique_ptr<TimeIntegrator> integrator = factory.createTimeIntegrator(Method, restricted3body_func, initial_value, step_size, order, num_iterations[i], Richardson);

            if (integrator == nullptr) 
            {
                std::cerr << "Error: Failed to create the specified time integrator." << std::endl;
                return 1;
            }

            std::cout << "---------" << Method << ", T = " << period_time << ", order = " << order << ", k = T/" << num_iterations[i] << ", with Richardson extrapolation" << "---------" << std::endl;

            auto start_time = std::chrono::high_resolution_clock::now();
            integrator->solve();
            auto end_time = std::chrono::high_resolution_clock::now();
            std::vector<std::vector<double>> results = integrator->get_Results();
            auto CPU_time = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time).count();
            std::vector<double> diff = integrator->get_error();
            std::vector<std::vector<double>> convergence_rate = integrator->get_convergence_rate();

            std::cout << "E[h]: ";
            for (const auto& value : convergence_rate[0])
            {
                std::cout << value << " ";
            }
            std::cout << std::endl;
            error_h_list.push_back(convergence_rate[0]);
            std::cout << "E[h/2]: ";
            for (const auto& value : convergence_rate[1])
            {
                std::cout << value << " ";
            }
            std::cout << std::endl;
            error_halfh_list.push_back(convergence_rate[1]);
            std::cout << "convergence rate: ";
            for (const auto& value : convergence_rate[2])
            {
                std::cout << value << " ";
            }
            std::cout << std::endl;
            convergence_rate_list.push_back(convergence_rate[2]);
            std::cout << "average convergence rate: ";
            std::cout << 0.25*(convergence_rate[2][0] + convergence_rate[2][1] + convergence_rate[2][3] + convergence_rate[2][4]) << std::endl;
            average_rate.push_back(0.25*(convergence_rate[2][0] + convergence_rate[2][1] + convergence_rate[2][3] + convergence_rate[2][4]));
            std::cout << "CPU time(ms): " << CPU_time << std::endl;
            time.push_back(CPU_time);

            if (output)
            {
                std::ofstream output_file("results.csv");
                output_file << "u1,u2,u3,u4,u5,u6\n";
                for (const auto& result : results) 
                {
                    for (size_t i = 0; i < result.size(); ++i) 
                    {
                        output_file << result[i];
                        if (i < result.size() - 1) 
                        {
                            output_file << ",";
                        }
                    }
                    output_file << "\n";
                }
                output_file.close();
            }
        }

        std::cout << "----------------------------LaTeX output----------------------------" << std::endl;
        std::vector<std::string> u_labels = {"u_1", "u_2", "u_3", "u_4", "u_5", "u_6"};
        for (size_t i = 0; i < u_labels.size(); ++i) 
        {
            std::cout << "$" << u_labels[i] << "$ ";
            for (size_t j = 0; j < error_h_list.size(); ++j) 
            {
                if (error_h_list[j][i] == 0.0)
                {
                    std::cout << "& " << 0 << "& " << 0 << " &$\\backslash$ " <<  " ";
                }
                else
                {
                    std::cout << "& " << std::scientific << std::setprecision(3) << error_h_list[j][i] << " &" << std::scientific << std::setprecision(3) << error_halfh_list[j][i] << " &" << std::fixed << std::setprecision(3) << convergence_rate_list[j][i] << " ";
                }     
            }
            std::cout << "\\\\" << std::endl;
        }   
        std::cout << "\\hline" << std::endl;
        std::cout << "平均收敛阶 ";
        for (size_t i = 0; i < average_rate.size(); ++i) 
        {
            if (i == average_rate.size()-1)
            {
                std::cout << "& \\multicolumn{3}{c}{" << std::fixed << std::setprecision(2) << average_rate[i] << "} ";
            }
            else
            {
                std::cout << "& \\multicolumn{3}{c|}{" << std::fixed << std::setprecision(2) << average_rate[i] << "} ";
            }
            
        }
        std::cout << "\\\\" << std::endl;

        std::cout << "\\hline" << std::endl;
        std::cout << "CPU time(ms) ";
        for (size_t i = 0; i < time.size(); ++i) 
        {
            if (i == time.size()-1)
            {
                std::cout << "& \\multicolumn{3}{c}{" <<  std::fixed << std::setprecision(0) << time[i] << "} ";
            }
            else
            {
                std::cout << "& \\multicolumn{3}{c|}{" <<  std::fixed << std::setprecision(0) << time[i] << "} ";
            }
        }
        std::cout << "\\\\" << std::endl;
        std::cout << "\\hline" << std::endl;
    }
    else
    {
        std::vector<std::vector<double>> error;
        std::vector<std::vector<double>> convergence_rate_list;
        std::vector<double> average_rate;
        std::vector<double> time;
        for (int i = 0; i < num_iterations.size(); ++i) 
        {
            TimeIntegratorFactory &factory = TimeIntegratorFactory::getInstance();
            double step_size = period_time / num_iterations[i];
            std::unique_ptr<TimeIntegrator> integrator = factory.createTimeIntegrator(Method, restricted3body_func, initial_value, step_size, order, num_iterations[i], Richardson);

            if (integrator == nullptr) 
            {
                std::cerr << "Error: Failed to create the specified time integrator." << std::endl;
                return 1;
            }

            std::cout << "--------------" << Method << ", T = " << period_time << ", order = " << order << ", k = T/" << num_iterations[i] << "--------------" << std::endl;

            auto start_time = std::chrono::high_resolution_clock::now();
            integrator->solve();
            auto end_time = std::chrono::high_resolution_clock::now();
            auto CPU_time = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time).count();

            std::vector<std::vector<double>> results = integrator->get_Results();
            std::vector<double> diff = integrator->get_error();

            std::cout << "absolute error :" << std::endl;
            for (const auto& value : diff)
            {
                std::cout << value << " ";
            }
            std::cout << std::endl;
            error.push_back(diff);

            if (output)
            {
                std::ofstream output_file("results.csv");
                output_file << "u1,u2,u3,u4,u5,u6\n";
                for (const auto& result : results) 
                {
                    for (size_t i = 0; i < result.size(); ++i) 
                    {
                        output_file << result[i];
                        if (i < result.size() - 1) 
                        {
                            output_file << ",";
                        }
                    }
                    output_file << "\n";
                }
                output_file.close();
            }

            std::vector<double> convergence_rate(error[0].size(), 0.0);
            if (i != 0)
            {
                for (size_t j = 0; j < error[i].size(); ++j)
                {
                    convergence_rate[j] = log2(error[i - 1][j] / error[i][j]);
                }
                convergence_rate_list.push_back(convergence_rate);
            }
            std::cout << "convergence rate: ";
            for (const auto& value : convergence_rate)
            {
                std::cout << value << " ";
            }
            std::cout << std::endl;
            std::cout << "average convergence rate: ";
            std::cout << 0.25*(convergence_rate[0] + convergence_rate[1] + convergence_rate[3] + convergence_rate[4]) << std::endl;
            average_rate.push_back(0.25*(convergence_rate[0] + convergence_rate[1] + convergence_rate[3] + convergence_rate[4]));
            std::cout << "CPU time(ms): " << CPU_time << std::endl;
            time.push_back(CPU_time);
        } 

        std::cout << "----------------------------LaTeX output----------------------------" << std::endl;
        std::vector<std::string> u_labels = {"u_1", "u_2", "u_3", "u_4", "u_5", "u_6"};
        for (size_t i = 0; i < u_labels.size(); ++i) 
        {
            std::cout << "$" << u_labels[i] << "$ ";
            for (size_t j = 0; j < error.size(); ++j) 
            {
                if (j == 0)
                {
                    if (error[j][i] == 0.0)
                    {
                        std::cout << "& " << 0 << " &$\\backslash$ " <<  " ";
                    }
                    else
                    {
                        std::cout << "& " << std::scientific << std::setprecision(3) << error[j][i] << " &$\\backslash$ " <<  " ";
                    }
                }
                else
                {
                    if (error[j][i] == 0.0)
                    {
                        std::cout << "& " << 0 << " &$\\backslash$ " <<  " ";
                    }
                    else
                    {
                       std::cout << "& " << std::scientific << std::setprecision(3) << error[j][i] << " &" << std::fixed << std::setprecision(3) << convergence_rate_list[j-1][i] << " ";
                    }
                    
                }     
            }
            std::cout << "\\\\" << std::endl;
        }   
        std::cout << "\\hline" << std::endl;
        std::cout << "平均收敛阶 ";
        std::cout << "& \\multicolumn{2}{c|}{" << " $\\backslash$ " << "} ";
        for (size_t i = 1; i < average_rate.size(); ++i) 
        {
            if (i == average_rate.size()-1)
            {
                std::cout << "& \\multicolumn{2}{c}{" << std::fixed << std::setprecision(2) << average_rate[i] << "} ";
            }
            else
            {
                std::cout << "& \\multicolumn{2}{c|}{" << std::fixed << std::setprecision(2) << average_rate[i] << "} ";
            }
            
        }
        std::cout << "\\\\" << std::endl;

        std::cout << "\\hline" << std::endl;
        std::cout << "CPU time(ms) ";
        for (size_t i = 0; i < time.size(); ++i) 
        {
            if (i == time.size()-1)
            {
                std::cout << "& \\multicolumn{2}{c}{" <<  std::fixed << std::setprecision(0) << time[i] << "} ";
            }
            else
            {
                std::cout << "& \\multicolumn{2}{c|}{" <<  std::fixed << std::setprecision(0) << time[i] << "} ";
            }
        }
        std::cout << "\\\\" << std::endl;
        std::cout << "\\hline" << std::endl;
    }

    return 0;
}



