#include <iostream>
#include <memory>
#include <string>
#include <map>
#include <functional>
#include <algorithm>
#include <vector>
#include <cmath>


class Vector_Function 
{
public:
    virtual std::vector<double> operator()(std::vector<double>& input, double time) 
    {
        return std::vector<double>(); 
    }
};

class TimeIntegrator 
{
public:
    virtual void solve() = 0;
    virtual std::vector<std::vector<double>> get_Results() = 0;
    virtual std::vector<double> get_error() = 0;
    virtual std::vector<std::vector<double>> get_convergence_rate() = 0;
};

class Adams_Bashforth : public TimeIntegrator 
{
public:
    Adams_Bashforth(Vector_Function& _vecFunc,const std::vector<double>& _initial_value, double _k, int _p, int _iteration_num, std::string _if_Richardson)
        : vecFunc(_vecFunc), initial_value(_initial_value), k(_k), p(_p), iteration_num(_iteration_num), if_Richardson(_if_Richardson)  
    {
        if (p < 1 || p > 4) 
        {
            throw std::invalid_argument("p must be between 1 and 4");
        }
        if (if_Richardson != "yes" && if_Richardson != "no") 
        {
            throw std::invalid_argument("invalid option for Richardson extrapolation");
        }

        // initialize beta
        switch (p) {
            case 1:
                betas = {1.0};
                break;
            case 2:
                betas = {-0.5, 1.5};
                break;
            case 3:
                betas = {5.0 / 12.0, -16.0 / 12.0, 23.0 / 12.0};
                break;
            case 4:
                betas = {-9.0 / 24.0, 37.0 / 24.0, -59.0 / 24.0, 55.0 / 24.0};
                break;
        }
    }

    void solve() override 
    {
        results = solve_Adams_Bashforth(k, iteration_num);
    }

    std::vector<std::vector<double>> get_Results() override
    {
        return results;
    }

    std::vector<std::vector<double>> get_convergence_rate() override
    {
        if (if_Richardson == "no")
        {
            throw std::invalid_argument("invalid use of get_convergence_rate without Richardson extrapolation");
        }

        std::vector<std::vector<double>> convergence_rate;
        results_half = solve_Adams_Bashforth(k * 0.5, iteration_num * 2);
        results_quat = solve_Adams_Bashforth(k * 0.25, iteration_num * 4);

        std::vector<double> diff_h(initial_value.size(), 0);
        std::vector<double> diff_halfh(initial_value.size(), 0);
        std::vector<double> rate(initial_value.size(), 0);
        for (size_t i = 0; i < initial_value.size(); ++i)
        {
            diff_h[i] = std::abs(results.back()[i] - results_quat.back()[i]);
            diff_halfh[i] = std::abs(results_half.back()[i] - results_quat.back()[i]);
            rate[i] = log2(diff_h[i] / diff_halfh[i] - 1.0);
        }

        convergence_rate.push_back(diff_h);
        convergence_rate.push_back(diff_halfh);
        convergence_rate.push_back(rate);

        return convergence_rate;
    }

    std::vector<double> get_error()
    {
        std::vector<double> diff(initial_value.size(), 0);
        if (!results.empty())
        {
            for (size_t i = 0; i < initial_value.size(); ++i)
            {
                diff[i] = std::abs(results.back()[i] - results.front()[i]);
            }
        }
        return diff;
    }

private:
    std::vector<std::vector<double>> solve_Adams_Bashforth(double curtime_step, int curiter_num)
    {
        std::vector<std::vector<double>> values(p, initial_value);
        std::vector<std::vector<double>> curresults;

        // Obtain sufficient initial values using 4th-order Runge-Kutta
        for (int i = 1; i < p; ++i) 
        {
            values[i] = runge_kutta_4_step(values[i - 1], curtime_step, i - 1);
        }

        // curresults stores the results of each iteration
        curresults.reserve(curiter_num + 1);

        for (const auto& value : values) 
        {
            curresults.push_back(value);
        }

        // Adam-Bashforth
        for (int i = 0; i < curiter_num - p + 1; ++i) 
        {
            std::vector<double> new_value(initial_value.size(), 0);
            std::vector<double> func_value;
            for (int j = 0; j < p; ++j) 
            {
                func_value = vecFunc(values[j], (i + j) * curtime_step);
                for (size_t m = 0; m < initial_value.size(); ++m) 
                {
                    new_value[m] += betas[j] * func_value[m];
                }
            }

            for (size_t m = 0; m < initial_value.size(); ++m) 
            {
                new_value[m] = values[p - 1][m] + curtime_step * new_value[m];
            }

            values.erase(values.begin());
            values.push_back(new_value);

            curresults.push_back(new_value);
        }

        return curresults;
    }

    std::vector<double> runge_kutta_4_step(std::vector<double>& value, double curtime_step, int n)
    {
        std::vector<double> k1 = vecFunc(value, n * curtime_step);
        std::vector<double> temp1(value.size());
        for (size_t i = 0; i < value.size(); ++i) 
        {
            temp1[i] = value[i] + curtime_step * k1[i] / 2;
        }

        std::vector<double> k2 = vecFunc(temp1, n * curtime_step + 0.5 * curtime_step);
        std::vector<double> temp2(value.size());
        for (size_t i = 0; i < value.size(); ++i) 
        {
            temp2[i] = value[i] + curtime_step * k2[i] / 2;
        }

        std::vector<double> k3 = vecFunc(temp2, n * curtime_step + 0.5 * curtime_step);
        std::vector<double> temp3(value.size());
        for (size_t i = 0; i < value.size(); ++i) 
        {
            temp3[i] = value[i] + curtime_step * k3[i];
        }

        std::vector<double> k4 = vecFunc(temp3, n * curtime_step + curtime_step);
        std::vector<double> result(value.size());
        for (size_t i = 0; i < value.size(); ++i) 
        {
            result[i] = value[i] + curtime_step * (k1[i] + 2 * k2[i] + 2 * k3[i] + k4[i]) / 6;
        }

        return result;
    }

    Vector_Function& vecFunc;
    const std::vector<double> initial_value;
    double k;
    int p;
    int iteration_num;
    std::vector<std::vector<double>> results;
    std::vector<std::vector<double>> results_half;
    std::vector<std::vector<double>> results_quat;
    std::string if_Richardson;
    std::vector<double> betas;
};

class Adams_Moulton : public TimeIntegrator
{
public:
    Adams_Moulton(Vector_Function& _vecFunc,const std::vector<double>& _initial_value, double _k, int _p, int _iteration_num, std::string _if_Richardson)
        : vecFunc(_vecFunc), initial_value(_initial_value), k(_k), p(_p), iteration_num(_iteration_num), if_Richardson(_if_Richardson)
    {
        if (p < 2 || p > 5)
        {
            throw std::invalid_argument("p must be between 2 and 5");
        }

        // initialize betas
        switch (p)
        {
            case 2:
                betas = {1.0 / 2.0, 1.0 / 2.0};
                break;
            case 3:
                betas = {-1.0 / 12.0, 8.0 / 12.0, 5.0 / 12.0};
                break;
            case 4:
                betas = {1.0 / 24.0, -5.0 / 24.0, 19.0 / 24.0, 9.0 / 24.0};
                break;
            case 5:
                betas = {-19.0 / 720.0, 106.0 / 720.0, -264.0 / 720.0, 646.0 / 720.0, 251.0 / 720.0};
                break;
        }
    }

    void solve() override
    {
        results = solve_Adams_Moulton(k, iteration_num);    
    }

    std::vector<std::vector<double>> get_Results() override
    {
        return results;
    }

    std::vector<std::vector<double>> get_convergence_rate() override
    {
        if (if_Richardson == "no")
        {
            throw std::invalid_argument("invalid use of get_convergence_rate without Richardson extrapolation");
        }

        std::vector<std::vector<double>> convergence_rate;
        results_half = solve_Adams_Moulton(k * 0.5, iteration_num * 2);
        results_quat = solve_Adams_Moulton(k * 0.25, iteration_num * 4);

        std::vector<double> diff_h(initial_value.size(), 0);
        std::vector<double> diff_halfh(initial_value.size(), 0);
        std::vector<double> rate(initial_value.size(), 0);
        for (size_t i = 0; i < initial_value.size(); ++i)
        {
            diff_h[i] = std::abs(results.back()[i] - results_quat.back()[i]);
            diff_halfh[i] = std::abs(results_half.back()[i] - results_quat.back()[i]);
            rate[i] = log2(diff_h[i] / diff_halfh[i] - 1.0);
        }

        convergence_rate.push_back(diff_h);
        convergence_rate.push_back(diff_halfh);
        convergence_rate.push_back(rate);

        return convergence_rate;
    }

    std::vector<double> get_error()
    {
        std::vector<double> diff(initial_value.size(), 0);
        if (!results.empty())
        {
            for (size_t i = 0; i < initial_value.size(); ++i)
            {
                diff[i] = std::abs(results.back()[i] - results.front()[i]);
            }
        }
        return diff;
    }

private:
    std::vector<std::vector<double>> solve_Adams_Moulton(double curtime_step, int curiter_num)
    {
        std::vector<std::vector<double>> values(p - 1, initial_value);
        std::vector<std::vector<double>> curresults;

        // Obtain sufficient initial values using 5th-order Runge-Kutta
        for (int i = 1; i < p - 1; ++i)
        {
            values[i] = runge_kutta_5_step(values[i - 1], curtime_step, i - 1);
        }

        // curresults stores the results of each iteration
        curresults.reserve(curiter_num + 1);

        for (const auto& value : values) 
        {
            curresults.push_back(value);
        }

        // Adam-Moulton
        for (int i = 0; i < curiter_num - p + 2; ++i)
        {
            std::vector<double> new_value(initial_value.size(), 0);
            new_value = fixed_point_iteration(values, i, curtime_step);

            values.erase(values.begin());
            values.push_back(new_value);

            curresults.push_back(new_value);
        }

        return curresults;
    }

    std::vector<double> runge_kutta_5_step(std::vector<double>& value, double curtime_step, int n)
    {
        std::vector<double> result(value.size());

        // Initialize RK matrix, RK weights, RK nodes
        std::vector<std::vector<double>> A = {
            {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
            {1.0 / 5.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
            {3.0 / 40.0, 9.0 / 40.0, 0.0, 0.0, 0.0, 0.0, 0.0},
            {44.0 / 45.0, -56.0 / 15.0, 32.0 / 9.0, 0.0, 0.0, 0.0, 0.0},
            {19372.0 / 6561.0, -25360.0 / 2187.0, 64448.0 / 6561.0, -212.0 / 729.0, 0.0, 0.0, 0.0},
            {9017.0 / 3168.0, -355.0 / 33.0, 46732.0 / 5247.0, 49.0 / 176.0, -5103.0 / 18656.0, 0.0, 0.0},
            {35.0 / 384.0, 0.0, 500.0 / 1113.0, 125.0 / 192.0, -2187.0 / 6784.0, 11.0 / 84.0, 0.0}
        };
        std::vector<double> b = {35.0 / 384.0, 0.0, 500.0 / 1113.0, 125.0 / 192.0, -2187.0 / 6784.0, 11.0 / 84.0, 0.0};
        std::vector<double> c = {0.0, 0.2, 0.3, 0.8, 8.0 / 9.0, 1.0, 1.0};

        int s = A.size();
        std::vector<std::vector<double>> y(s, std::vector<double>(value.size(), 0));

        for (int j = 0; j < s; ++j)
        {
            std::vector<double> U_temp = value;
            for (int l = 0; l < j; ++l)
            {
                for (size_t m = 0; m < U_temp.size(); ++m)
                {
                    U_temp[m] += curtime_step * A[j][l] * y[l][m];
                }
            }
            y[j] = vecFunc(U_temp, n * curtime_step + c[j] * curtime_step);
        }

        result = value;
        for (int j= 0; j < s; ++j)
        {
            for (size_t m = 0; m < value.size(); ++m)
            {
                result[m] += curtime_step * b[j] * y[j][m];
            }
        }
        return result;
    }

    std::vector<double> fixed_point_iteration(std::vector<std::vector<double>>& values, int n, double time_step)
    {
        const double tol = 1e-15;
        const int max_iter = 100;
        std::vector<double> prev_value = values.back();
        bool converged;

        for (int iter = 0; iter < max_iter; ++iter)
        {
            converged = true;
            std::vector<double> curr_value(prev_value.size(), 0);
            std::vector<double> func_value;
            for (int j = 0; j < p - 1; ++j) 
            {
                func_value = vecFunc(values[j], (n + j) * time_step);
                for (size_t m = 0; m < curr_value.size(); ++m) 
                {
                    curr_value[m] += betas[j] * func_value[m];
                }
            }

            std::vector<double> pre_func = vecFunc(prev_value, (n + p - 2) * time_step);
            for (size_t m = 0; m < curr_value.size(); ++m) 
            {
                curr_value[m] += betas[p - 1] * pre_func[m];
            }

            for (size_t m = 0; m < curr_value.size(); ++m) 
            {
                curr_value[m] = values[p - 2][m] + time_step * curr_value[m];
            }

            for (size_t m = 0; m < curr_value.size(); ++m)
            {
                if(std::abs(curr_value[m] - prev_value[m]) > tol)
                {
                    converged = false;
                }
            }

            if (converged)
            {
                return curr_value;
            }
            prev_value = curr_value;
        }

        std::cerr << "Warning: fixed_point_iteration did not converge within the maximum number of iterations." << std::endl;
        return prev_value;
    }

    Vector_Function& vecFunc;
    std::vector<double> initial_value;
    double k;
    int p;
    int iteration_num;
    std::vector<std::vector<double>> results;
    std::vector<std::vector<double>> results_half;
    std::vector<std::vector<double>> results_quat;
    std::string if_Richardson;
    std::vector<double> betas;
};

class Backward_differentiation : public TimeIntegrator
{
public:
    Backward_differentiation(Vector_Function& _vecFunc,const std::vector<double>& _initial_value, double _k, int _p, int _iteration_num, std::string _if_Richardson)
        : vecFunc(_vecFunc), initial_value(_initial_value), k(_k), p(_p), iteration_num(_iteration_num), if_Richardson(_if_Richardson)
    {
        if (p < 1 || p > 4)
        {
            throw std::invalid_argument("p must be between 2 and 5");
        }

        // initialize alphas, beta
        switch (p)
        {  
            case 1:
                alphas = {-1.0};
                beta = 1.0;
                break;
            case 2:
                alphas = {1.0 / 3.0, -4.0 / 3.0};
                beta = 2.0 / 3.0;
                break;
            case 3:
                alphas = {-2.0 / 11.0, 9.0 / 11.0, -18.0 / 11.0};
                beta = 6.0 / 11.0;
                break;
            case 4:
                alphas = {3.0 / 25.0, -16.0 / 25.0, 36.0 / 25.0, -48.0 / 25.0};
                beta = 12.0 / 25.0;
                break;
        }
    }

    void solve() override
    {
        results = solve_Backward_differentiation(k, iteration_num);   
    }

    std::vector<std::vector<double>> get_Results() override
    {
        return results;
    }

    std::vector<std::vector<double>> get_convergence_rate() override
    {
        if (if_Richardson == "no")
        {
            throw std::invalid_argument("invalid use of get_convergence_rate without Richardson extrapolation");
        }

        std::vector<std::vector<double>> convergence_rate;
        results_half = solve_Backward_differentiation(k * 0.5, iteration_num * 2);
        results_quat = solve_Backward_differentiation(k * 0.25, iteration_num * 4);

        std::vector<double> diff_h(initial_value.size(), 0);
        std::vector<double> diff_halfh(initial_value.size(), 0);
        std::vector<double> rate(initial_value.size(), 0);
        for (size_t i = 0; i < initial_value.size(); ++i)
        {
            diff_h[i] = std::abs(results.back()[i] - results_quat.back()[i]);
            diff_halfh[i] = std::abs(results_half.back()[i] - results_quat.back()[i]);
            rate[i] = log2(diff_h[i] / diff_halfh[i] - 1.0);
        }

        convergence_rate.push_back(diff_h);
        convergence_rate.push_back(diff_halfh);
        convergence_rate.push_back(rate);

        return convergence_rate;
    }

    std::vector<double> get_error()
    {
        std::vector<double> diff(initial_value.size(), 0);
        if (!results.empty())
        {
            for (size_t i = 0; i < initial_value.size(); ++i)
            {
                diff[i] = std::abs(results.back()[i] - results.front()[i]);
            }
        }
        return diff;
    }

private:
    std::vector<std::vector<double>> solve_Backward_differentiation(double curtime_step, int curiter_num)
    {
        std::vector<std::vector<double>> values(p, initial_value);
        std::vector<std::vector<double>> curresults;

        // Obtain sufficient initial values using 4th-order Runge-Kutta
        for (int i = 1; i < p; ++i)
        {
            values[i] = runge_kutta_4_step(values[i - 1], curtime_step, i - 1);
        }

        // curresults stores the results of each iteration
        curresults.reserve(curiter_num + 1);

        for (const auto& value : values)
        {
            curresults.push_back(value);
        }

        // BDF
        for (int i = 0; i < curiter_num - p + 1; ++i)
        {
            std::vector<double> new_value(initial_value.size(), 0);
            new_value = fixed_point_iteration(values, i, curtime_step);

            values.erase(values.begin());
            values.push_back(new_value);

            curresults.push_back(new_value);
        }

        return curresults;
    }

    std::vector<double> runge_kutta_4_step(std::vector<double>& value, double curtime_step, int n)
    {
        std::vector<double> k1 = vecFunc(value, n * curtime_step);
        std::vector<double> temp1(value.size());
        for (size_t i = 0; i < value.size(); ++i) 
        {
            temp1[i] = value[i] + curtime_step * k1[i] / 2;
        }

        std::vector<double> k2 = vecFunc(temp1, n * curtime_step + 0.5 * curtime_step);
        std::vector<double> temp2(value.size());
        for (size_t i = 0; i < value.size(); ++i) 
        {
            temp2[i] = value[i] + curtime_step * k2[i] / 2;
        }

        std::vector<double> k3 = vecFunc(temp2, n * curtime_step + 0.5 * curtime_step);
        std::vector<double> temp3(value.size());
        for (size_t i = 0; i < value.size(); ++i) 
        {
            temp3[i] = value[i] + curtime_step * k3[i];
        }

        std::vector<double> k4 = vecFunc(temp3, n * curtime_step + curtime_step);
        std::vector<double> result(value.size());
        for (size_t i = 0; i < value.size(); ++i) 
        {
            result[i] = value[i] + curtime_step * (k1[i] + 2 * k2[i] + 2 * k3[i] + k4[i]) / 6;
        }

        return result;
    }

    std::vector<double> fixed_point_iteration(std::vector<std::vector<double>>& values, int n, double time_step)
    {
        const double tol = 1e-15;
        const int max_iter = 100;
        std::vector<double> prev_value = values.back();
        bool converged;

        for (int iter = 0; iter < max_iter; ++iter)
        {
            converged = true;
            std::vector<double> curr_value(prev_value.size(), 0);
            for (int j = 0; j < p; ++j) 
            {
                for (size_t m = 0; m < curr_value.size(); ++m) 
                {
                    curr_value[m] += -alphas[j] * values[j][m];
                }
            }

            std::vector<double> pre_func = vecFunc(prev_value, (n + p - 1) * time_step);
            for (size_t m = 0; m < curr_value.size(); ++m) 
            {
                curr_value[m] += time_step * beta * pre_func[m];
            }

            for (size_t m = 0; m < curr_value.size(); ++m)
            {
                if(std::abs(curr_value[m] - prev_value[m]) > tol)
                {
                    converged = false;
                }
            }

            if (converged)
            {
                return curr_value;
            }
            prev_value = curr_value;
        }

        std::cerr << "Warning: fixed_point_iteration did not converge within the maximum number of iterations." << std::endl;
        return prev_value;
    }

    Vector_Function& vecFunc;
    std::vector<double> initial_value;
    double k;
    int p;
    int iteration_num;
    std::vector<std::vector<double>> results;
    std::vector<std::vector<double>> results_half;
    std::vector<std::vector<double>> results_quat;
    std::string if_Richardson;
    std::vector<double> alphas;
    double beta;
};

class Classical_RK : public TimeIntegrator
{
public:
    public:
    Classical_RK(Vector_Function& _vecFunc,const std::vector<double>& _initial_value, double _k, int _p, int _iteration_num, std::string _if_Richardson)
        : vecFunc(_vecFunc), initial_value(_initial_value), k(_k), p(_p), iteration_num(_iteration_num), if_Richardson(_if_Richardson)
    {
        if (p != 4)
        {
            throw std::invalid_argument("p must be 4");
        }
    }

    void solve() override
    {
        results = solve_Classical_RK(k, iteration_num); 
    }

    std::vector<std::vector<double>> get_Results() override
    {
        return results;
    }

    std::vector<std::vector<double>> get_convergence_rate() override
    {
        if (if_Richardson == "no")
        {
            throw std::invalid_argument("invalid use of get_convergence_rate without Richardson extrapolation");
        }

        std::vector<std::vector<double>> convergence_rate;
        results_half = solve_Classical_RK(k * 0.5, iteration_num * 2);
        results_quat = solve_Classical_RK(k * 0.25, iteration_num * 4);

        std::vector<double> diff_h(initial_value.size(), 0);
        std::vector<double> diff_halfh(initial_value.size(), 0);
        std::vector<double> rate(initial_value.size(), 0);
        for (size_t i = 0; i < initial_value.size(); ++i)
        {
            diff_h[i] = std::abs(results.back()[i] - results_quat.back()[i]);
            diff_halfh[i] = std::abs(results_half.back()[i] - results_quat.back()[i]);
            rate[i] = log2(diff_h[i] / diff_halfh[i] - 1.0);
        }

        convergence_rate.push_back(diff_h);
        convergence_rate.push_back(diff_halfh);
        convergence_rate.push_back(rate);

        return convergence_rate;
    }

    std::vector<double> get_error()
    {
        std::vector<double> diff(initial_value.size(), 0);
        if (!results.empty())
        {
            for (size_t i = 0; i < initial_value.size(); ++i)
            {
                diff[i] = std::abs(results.back()[i] - results.front()[i]);
            }
        }
        return diff;
    }

private:
    std::vector<std::vector<double>> solve_Classical_RK(double curtime_step, int curiter_num)
    {
        std::vector<std::vector<double>> curresults;
        curresults.push_back(initial_value);

        for (int i = 0; i < curiter_num; ++i)
        {
            std::vector<double> new_value = runge_kutta_4_step(curresults.back(), curtime_step, i);
            curresults.push_back(new_value);
        }

        return curresults;
    }

    std::vector<double> runge_kutta_4_step(std::vector<double>& value, double curtime_step, int n)
    {
        std::vector<double> k1 = vecFunc(value, n * curtime_step);
        std::vector<double> temp1(value.size());
        for (size_t i = 0; i < value.size(); ++i) 
        {
            temp1[i] = value[i] + curtime_step * k1[i] / 2;
        }

        std::vector<double> k2 = vecFunc(temp1, n * curtime_step + 0.5 * curtime_step);
        std::vector<double> temp2(value.size());
        for (size_t i = 0; i < value.size(); ++i) 
        {
            temp2[i] = value[i] + curtime_step * k2[i] / 2;
        }

        std::vector<double> k3 = vecFunc(temp2, n * curtime_step + 0.5 * curtime_step);
        std::vector<double> temp3(value.size());
        for (size_t i = 0; i < value.size(); ++i) 
        {
            temp3[i] = value[i] + curtime_step * k3[i];
        }

        std::vector<double> k4 = vecFunc(temp3, n * curtime_step + curtime_step);
        std::vector<double> result(value.size());
        for (size_t i = 0; i < value.size(); ++i) 
        {
            result[i] = value[i] + curtime_step * (k1[i] + 2 * k2[i] + 2 * k3[i] + k4[i]) / 6;
        }

        return result;
    }

    Vector_Function& vecFunc;
    std::vector<double> initial_value;
    double k;
    int p;
    int iteration_num;
    std::vector<std::vector<double>> results;
    std::vector<std::vector<double>> results_half;
    std::vector<std::vector<double>> results_quat;
    std::string if_Richardson;
};

class explicit_SDIRK : public TimeIntegrator
{
public:
    public:
    explicit_SDIRK(Vector_Function& _vecFunc,const std::vector<double>& _initial_value, double _k, int _p, int _iteration_num, std::string _if_Richardson)
        : vecFunc(_vecFunc), initial_value(_initial_value), k(_k), p(_p), iteration_num(_iteration_num), if_Richardson(_if_Richardson)
    {
        if (p != 4)
        {
            throw std::invalid_argument("p must be 4");
        }

        A = {
            {0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
            {1.0 / 4.0, 1.0 / 4.0, 0.0, 0.0, 0.0, 0.0},
            {8611.0 / 62500.0, -1743.0 / 31250.0, 1.0 / 4.0, 0.0, 0.0, 0.0},
            {5012029.0 / 34652500.0, -654441.0 / 2922500.0, 174375.0 / 388108.0, 1.0 / 4.0, 0.0, 0.0},
            {15267082809.0 / 155376265600.0, -71443401.0 / 120774400.0, 730878875.0 / 902184768.0, 2285395.0 / 8070912.0, 1.0 / 4.0, 0.0},
            {82889.0 / 524892.0, 0.0, 15625.0 / 83664.0, 69875.0 / 102672.0, -2260.0 / 8211.0, 1.0 / 4.0}
        };

        b = {82889.0 / 524892.0, 0.0, 15625.0 / 83664.0, 69875.0 / 102672.0, -2260.0 / 8211.0, 1.0 / 4.0};

        c = {0.0, 1.0 / 2.0, 83.0 / 250.0, 31.0 / 50.0, 17.0 / 20.0, 1.0};
    }

    void solve() override
    {
        results = solve_explicit_SDIRK(k, iteration_num);
    }

    std::vector<std::vector<double>> get_Results() override
    {
        return results;
    }

    std::vector<std::vector<double>> get_convergence_rate() override
    {
        if (if_Richardson == "no")
        {
            throw std::invalid_argument("invalid use of get_convergence_rate without Richardson extrapolation");
        }

        std::vector<std::vector<double>> convergence_rate;
        results_half = solve_explicit_SDIRK(k * 0.5, iteration_num * 2);
        results_quat = solve_explicit_SDIRK(k * 0.25, iteration_num * 4);

        std::vector<double> diff_h(initial_value.size(), 0);
        std::vector<double> diff_halfh(initial_value.size(), 0);
        std::vector<double> rate(initial_value.size(), 0);
        for (size_t i = 0; i < initial_value.size(); ++i)
        {
            diff_h[i] = std::abs(results.back()[i] - results_quat.back()[i]);
            diff_halfh[i] = std::abs(results_half.back()[i] - results_quat.back()[i]);
            rate[i] = log2(diff_h[i] / diff_halfh[i] - 1.0);
        }

        convergence_rate.push_back(diff_h);
        convergence_rate.push_back(diff_halfh);
        convergence_rate.push_back(rate);

        return convergence_rate;
    }

    std::vector<double> get_error()
    {
        std::vector<double> diff(initial_value.size(), 0);
        if (!results.empty())
        {
            for (size_t i = 0; i < initial_value.size(); ++i)
            {
                diff[i] = std::abs(results.back()[i] - results.front()[i]);
            }
        }
        return diff;
    }


private:
    std::vector<std::vector<double>> solve_explicit_SDIRK(double curtime_step, int curiter_num)
    {
        std::vector<std::vector<double>> curresults;
        curresults.push_back(initial_value);

        for (int i = 0; i < curiter_num; ++i)
        {
            std::vector<double> new_value(initial_value.size(), 0);
            new_value = fixed_point_iteration(curresults.back(), i, curtime_step);
            curresults.push_back(new_value);
        }

        return curresults;
    }

    std::vector<double> fixed_point_iteration(std::vector<double>& value, int n, double time_step)
    {
        const double tol = 1e-15;
        const int max_iter = 100;
        int s = A.size();
        std::vector<double> f_U_tn = vecFunc(value, n * time_step);
        std::vector<std::vector<double>> y(s, f_U_tn); 
        bool converged;

        for (int i = 0; i < s; ++i)
        {
            converged = false;
            for (int iter = 0; iter < max_iter; ++iter)
            {
                std::vector<double> U_temp = value;
                for (int j = 0; j <= i; ++j)
                {
                    for (int m = 0; m < value.size(); ++m)
                    {
                        U_temp[m] += time_step * A[i][j] * y[j][m];
                    }
                }

                std::vector<double> y_new = vecFunc(U_temp, n * time_step + c[i] * time_step);

                double error = 0;
                for (int j = 0; j < value.size(); ++j)
                {
                    error += std::abs(y_new[j] - y[i][j]);
                }

                if (error < tol)
                {
                    y[i] = y_new;
                    converged = true;
                    break;
                }

                y[i] = y_new;
            }
            if (!converged)
            {
                std::cerr << "Warning: fixed_point_iteration did not converge within the maximum number of iterations." << std::endl;
            }
        }

        std::vector<double> new_value(value.size(), 0);
        for (int i = 0; i < s; ++i)
        {
            for (int j = 0; j < value.size(); ++j)
            {
                new_value[j] += time_step * b[i] * y[i][j];
            }
        }

        for (int j = 0; j < value.size(); ++j)
        {
            new_value[j] += value[j];
        }

        return new_value;
    }

    Vector_Function& vecFunc;
    std::vector<double> initial_value;
    double k;
    int p;
    int iteration_num;
    std::vector<std::vector<double>> results;
    std::vector<std::vector<double>> results_half;
    std::vector<std::vector<double>> results_quat;
    std::string if_Richardson;
    std::vector<std::vector<double>> A;
    std::vector<double> b;
    std::vector<double> c;
};

class Gauss_Legendre_RK : public TimeIntegrator
{
public:
    public:
    Gauss_Legendre_RK(Vector_Function& _vecFunc,const std::vector<double>& _initial_value, double _k, int _p, int _iteration_num, std::string _if_Richardson)
        : vecFunc(_vecFunc), initial_value(_initial_value), k(_k), p(_p), iteration_num(_iteration_num), if_Richardson(_if_Richardson)
    {
        if (p != 2 && p != 4 && p != 6)
        {
            throw std::invalid_argument("p must be 2, 4 or 6");
        }

        // Initialize RK matrix, RK weights, RK nodes
        switch (p)
        {  
            case 2:
                A = {
                    {1.0 / 2.0}
                };
                b = {1.0};
                c = {1.0 / 2.0};
                break;
            case 4:
                A = {
                    {1.0 / 4.0, (3.0 - 2.0 * sqrt(3.0)) / 12.0},
                    {(3.0 + 2.0 * sqrt(3.0)) / 12.0, 1.0 / 4.0}
                };
                b = {1.0 / 2.0, 1.0 / 2.0};
                c = {(3.0 - sqrt(3.0)) / 6.0, (3.0 + sqrt(3.0)) / 6.0};
                break;
            case 6:
                A = {
                    {5.0 / 36.0, 2.0 / 9.0 - sqrt(15.0) / 15.0, 5.0 / 36.0 - sqrt(15.0) / 30.0},
                    {5.0 / 36.0 + sqrt(15.0) / 24.0, 2.0 / 9.0, 5.0 / 36.0 - sqrt(15.0) / 24.0},
                    {5.0 / 36.0 + sqrt(15.0) / 30.0, 2.0 / 9.0 + sqrt(15.0) / 15.0, 5.0 / 36.0}
                };
                b = {5.0 / 18.0, 4.0 / 9.0, 5.0 / 18.0};
                c = {(5.0 - sqrt(15.0)) / 10.0, 1.0 / 2.0, (5.0 + sqrt(15.0)) / 10.0};
                break;
        }
    }

    void solve() override
    {
        results = solve_Gauss_Legendre_RK(k, iteration_num);
    }

    std::vector<std::vector<double>> get_Results() override
    {
        return results;
    }

    std::vector<std::vector<double>> get_convergence_rate() override
    {
        if (if_Richardson == "no")
        {
            throw std::invalid_argument("invalid use of get_convergence_rate without Richardson extrapolation");
        }

        std::vector<std::vector<double>> convergence_rate;
        results_half = solve_Gauss_Legendre_RK(k * 0.5, iteration_num * 2);
        results_quat = solve_Gauss_Legendre_RK(k * 0.25, iteration_num * 4);

        std::vector<double> diff_h(initial_value.size(), 0);
        std::vector<double> diff_halfh(initial_value.size(), 0);
        std::vector<double> rate(initial_value.size(), 0);
        for (size_t i = 0; i < initial_value.size(); ++i)
        {
            diff_h[i] = std::abs(results.back()[i] - results_quat.back()[i]);
            diff_halfh[i] = std::abs(results_half.back()[i] - results_quat.back()[i]);
            rate[i] = log2(diff_h[i] / diff_halfh[i] - 1.0);
        }

        convergence_rate.push_back(diff_h);
        convergence_rate.push_back(diff_halfh);
        convergence_rate.push_back(rate);

        return convergence_rate;
    }

    std::vector<double> get_error()
    {
        std::vector<double> diff(initial_value.size(), 0);
        if (!results.empty())
        {
            for (size_t i = 0; i < initial_value.size(); ++i)
            {
                diff[i] = std::abs(results.back()[i] - results.front()[i]);
            }
        }
        return diff;
    }


private:
    std::vector<std::vector<double>> solve_Gauss_Legendre_RK(double curtime_step, int curiter_num)
    {
        std::vector<std::vector<double>> curresults;
        curresults.push_back(initial_value);

        for (int i = 0; i < curiter_num; ++i)
        {
            std::vector<double> new_value(initial_value.size(), 0);
            new_value = fixed_point_iteration(curresults.back(), i, curtime_step);
            curresults.push_back(new_value);
        }

        return curresults;
    }

    std::vector<double> fixed_point_iteration(std::vector<double>& value, int n, double time_step)
    {
        const double tol = 1e-15;
        const int max_iter = 100;
        int s = A.size();
        std::vector<double> f_U_tn = vecFunc(value, n * time_step);
        std::vector<std::vector<double>> y(s, f_U_tn); 
        bool converged;

        for (int iter = 0; iter < max_iter; ++iter)
        {
            std::vector<std::vector<double>> y_new(s, std::vector<double>(f_U_tn.size(), 0));
            converged = true;

            for (int i = 0; i < s; ++i)
            {
                std::vector<double> U_temp = value;
                for (int j = 0; j < s; ++j)
                {
                    for (size_t m = 0; m < value.size(); ++m)
                    {
                        U_temp[m] += time_step * A[i][j] * y[j][m];
                    }
                }

                y_new[i] = vecFunc(U_temp, n * time_step + c[i] * time_step);
            }

            for (int i = 0; i < s; ++i)
            {
                double error = 0;
                for (int j = 0; j < value.size(); ++j)
                {
                    error += std::abs(y_new[i][j] - y[i][j]);
                }

                if (error > tol)
                {
                    converged = false;
                }
            }
            
            if (converged == true)
            {
                y = y_new;
                break;
            }

            y = y_new;
        }

        if (!converged)
        {
            std::cerr << "Warning: fixed_point_iteration did not converge within the maximum number of iterations." << std::endl;
        }

        std::vector<double> new_value(value.size(), 0);
        for (int i = 0; i < s; ++i)
        {
            for (int j = 0; j < value.size(); ++j)
            {
                new_value[j] += time_step * b[i] * y[i][j];
            }
        }

        for (int j = 0; j < value.size(); ++j)
        {
            new_value[j] += value[j];
        }

        return new_value;
    }

    Vector_Function& vecFunc;
    std::vector<double> initial_value;
    double k;
    int p;
    int iteration_num;
    std::vector<std::vector<double>> results;
    std::vector<std::vector<double>> results_half;
    std::vector<std::vector<double>> results_quat;
    std::string if_Richardson;
    std::vector<std::vector<double>> A;
    std::vector<double> b;
    std::vector<double> c;
};

class Fehlberg_embedded_RK : public TimeIntegrator
{
public:
    public:
    Fehlberg_embedded_RK(Vector_Function& _vecFunc,const std::vector<double>& _initial_value, double _k, int _p, int _iteration_num, std::string _if_Richardson)
        : vecFunc(_vecFunc), initial_value(_initial_value), k(_k), p(_p), iteration_num(_iteration_num), if_Richardson(_if_Richardson)
    {
        if (p != 4)
        {
            throw std::invalid_argument("p must be 4 to ensure that the method is 4(5) embedded");
        }

        // Initialize RK matrix, RK weights, RK nodes
        A = {
            {0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
            {1.0 / 4.0, 0.0, 0.0, 0.0, 0.0, 0.0},
            {3.0 / 32.0, 9.0 / 32.0, 0.0, 0.0, 0.0, 0.0},
            {1932.0 / 2197.0, -7200.0 / 2197.0, 7296.0 / 2197.0, 0.0, 0.0, 0.0},
            {439.0 / 216.0, -8.0, 3680.0 / 513.0, -845.0 / 4104.0, 0.0, 0.0},
            {-8.0 / 27.0, 2.0, -3544.0 / 2565.0, 1859.0 / 4104.0, -11.0 / 40.0, 0.0}
        };
        b = {25.0 / 216.0, 0.0, 1408.0 / 2565.0, 2197.0 / 4104.0, -0.2, 0.0};
        b_hat = {16.0 / 135.0, 0.0, 6656.0 / 12825.0, 28561.0 / 56430.0, -9.0 / 50.0, 2.0 / 55.0};
        c = {0, 0.25, 0.375, 12.0 / 13.0, 1.0, 0.5};
    }

    void solve() override
    {
        int s = A.size();

        results.push_back(initial_value);
        double k_init = k;
        double T = 0.0;

        for (int i = 0; i < iteration_num; ++i)
        {
            std::vector<double> U_n = results.back();
            std::vector<double> U_np1, Uhat_np1;
            bool accepted = false;

            while (!accepted)
            {
                std::vector<std::vector<double>> y(s, std::vector<double>(U_n.size(), 0));

                for (int j = 0; j < s; ++j)
                {
                    std::vector<double> U_temp = U_n;
                    for (int l = 0; l < j; ++l)
                    {
                        for (size_t m = 0; m < U_temp.size(); ++m)
                        {
                            U_temp[m] += k_init * A[j][l] * y[l][m];
                        }
                    }
                    y[j] = vecFunc(U_temp, i * k_init + c[j] * k_init);
                }

                // Compute U^{n+1} and \hat{U}^{n+1}
                U_np1 = U_n;
                Uhat_np1 = U_n;
                for (int j= 0; j < s; ++j)
                {
                    for (size_t m = 0; m < U_n.size(); ++m)
                    {
                        U_np1[m] += k_init * b[j] * y[j][m];
                        Uhat_np1[m] += k_init * b_hat[j] * y[j][m];
                    }
                }

                // Compute error indicator E_ind
                double E_ind = 0.0;
                for (size_t m = 0; m < U_n.size(); ++m)
                {
                    double epsilon_m = E_abs + std::abs(U_n[m]) * E_rel;
                    E_ind += std::pow((Uhat_np1[m] - U_np1[m]) / epsilon_m, 2);
                }
                E_ind = std::sqrt(E_ind / U_n.size());

                // Compute new time step k_new
                double k_new = k_init * std::min(rho_max, std::max(rho_min, rho * std::pow(E_ind, -1.0 / 4.0)));

                if (T + k_init == k * iteration_num)
                {
                    results.push_back(U_np1);
                    break;
                }
                // Determine whether to accept the current step
                if (E_ind <= 1.0)
                {
                    if (T + k_new > k * iteration_num)
                    {
                        k_new = k * iteration_num - T;
                    }
                    else
                    {
                        accepted = true;
                        results.push_back(U_np1);
                    }
                }
                k_init = k_new;
            }

            T = T + k_init;
            if (T + k_init == k * iteration_num)
            {
                std::cout << "real step number is " << i << std::endl;
                break;
            }
        }
    }

    std::vector<std::vector<double>> get_Results() override
    {
        return results;
    }

    std::vector<std::vector<double>> get_convergence_rate() override
    {
        throw std::invalid_argument("convergence rate need not to be calculated in embedded RK");

        return results;
    }

    std::vector<double> get_error()
    {
        std::vector<double> diff(initial_value.size(), 0);
        if (!results.empty())
        {
            for (size_t i = 0; i < initial_value.size(); ++i)
            {
                diff[i] = std::abs(results.back()[i] - results.front()[i]);
            }
        }
        return diff;
    }


private:
    const double E_abs = 1e-8; 
    const double E_rel = 1e-8; 
    const double rho_max = 3;
    const double rho = 0.8;
    const double rho_min = 0.5;
    Vector_Function& vecFunc;
    std::vector<double> initial_value;
    double k;
    int p;
    int iteration_num;
    std::vector<std::vector<double>> results;
    std::string if_Richardson;
    std::vector<std::vector<double>> A;
    std::vector<double> b;
    std::vector<double> b_hat;
    std::vector<double> c;
};

class Dormand_Prince_embedded_RK : public TimeIntegrator
{
public:
    public:
    Dormand_Prince_embedded_RK(Vector_Function& _vecFunc,const std::vector<double>& _initial_value, double _k, int _p, int _iteration_num, std::string _if_Richardson)
        : vecFunc(_vecFunc), initial_value(_initial_value), k(_k), p(_p), iteration_num(_iteration_num), if_Richardson(_if_Richardson)
    {
        if (p != 5)
        {
            throw std::invalid_argument("p must be 5 to ensure that the method is 5(4) embedded");
        }

        // Initialize RK matrix, RK weights, RK nodes
        A = {
            {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
            {1.0 / 5.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
            {3.0 / 40.0, 9.0 / 40.0, 0.0, 0.0, 0.0, 0.0, 0.0},
            {44.0 / 45.0, -56.0 / 15.0, 32.0 / 9.0, 0.0, 0.0, 0.0, 0.0},
            {19372.0 / 6561.0, -25360.0 / 2187.0, 64448.0 / 6561.0, -212.0 / 729.0, 0.0, 0.0, 0.0},
            {9017.0 / 3168.0, -355.0 / 33.0, 46732.0 / 5247.0, 49.0 / 176.0, -5103.0 / 18656.0, 0.0, 0.0},
            {35.0 / 384.0, 0.0, 500.0 / 1113.0, 125.0 / 192.0, -2187.0 / 6784.0, 11.0 / 84.0, 0.0}
        };
        b = {35.0 / 384.0, 0.0, 500.0 / 1113.0, 125.0 / 192.0, -2187.0 / 6784.0, 11.0 / 84.0, 0.0};
        b_hat = {5179.0 / 57600.0, 0.0, 7571.0 / 16695.0, 393.0 / 640.0, -92097.0 / 339200.0, 187.0 / 2100.0, 1.0 / 40.0};
        c = {0.0, 0.2, 0.3, 0.8, 8.0 / 9.0, 1.0, 1.0};
    }

    void solve() override
    {
        int s = A.size();

        results.push_back(initial_value);
        double k_init = k;
        double T = 0.0;

        for (int i = 0; i < iteration_num; ++i)
        {
            std::vector<double> U_n = results.back();
            std::vector<double> U_np1, Uhat_np1;
            bool accepted = false;

            while (!accepted)
            {
                std::vector<std::vector<double>> y(s, std::vector<double>(U_n.size(), 0));

                for (int j = 0; j < s; ++j)
                {
                    std::vector<double> U_temp = U_n;
                    for (int l = 0; l < j; ++l)
                    {
                        for (size_t m = 0; m < U_temp.size(); ++m)
                        {
                            U_temp[m] += k_init * A[j][l] * y[l][m];
                        }
                    }
                    y[j] = vecFunc(U_temp, i * k_init + c[j] * k_init);
                }

                // Compute U^{n+1} and \hat{U}^{n+1}
                U_np1 = U_n;
                Uhat_np1 = U_n;
                for (int j= 0; j < s; ++j)
                {
                    for (size_t m = 0; m < U_n.size(); ++m)
                    {
                        U_np1[m] += k_init * b[j] * y[j][m];
                        Uhat_np1[m] += k_init * b_hat[j] * y[j][m];
                    }
                }

                // Compute error indicator E_ind
                double E_ind = 0.0;
                for (size_t m = 0; m < U_n.size(); ++m)
                {
                    double epsilon_m = E_abs + std::abs(U_n[m]) * E_rel;
                    E_ind += std::pow((Uhat_np1[m] - U_np1[m]) / epsilon_m, 2);
                }
                E_ind = std::sqrt(E_ind / U_n.size());

                // Compute new time step k_new
                double k_new = k_init * std::min(rho_max, std::max(rho_min, rho * std::pow(E_ind, -1.0 / 4.0)));

                if (T + k_init == k * iteration_num)
                {
                    results.push_back(U_np1);
                    break;
                }
                // Determine whether to accept the current step
                if (E_ind <= 1.0)
                {
                    if (T + k_new > k * iteration_num)
                    {
                        k_new = k * iteration_num - T;
                    }
                    else
                    {
                        accepted = true;
                        results.push_back(U_np1);
                    }
                }
                k_init = k_new;
            }

            T = T + k_init;
            if (T + k_init == k * iteration_num)
            {
                std::cout << "real step number is " << i << std::endl;
                break;
            }
        }
    }

    std::vector<std::vector<double>> get_Results() override
    {
        return results;
    }

    std::vector<std::vector<double>> get_convergence_rate() override
    {
        throw std::invalid_argument("convergence rate need not to be calculated in embedded RK");
        return results;
    }

    std::vector<double> get_error()
    {
        std::vector<double> diff(initial_value.size(), 0);
        if (!results.empty())
        {
            for (size_t i = 0; i < initial_value.size(); ++i)
            {
                diff[i] = std::abs(results.back()[i] - results.front()[i]);
            }
        }
        return diff;
    }


private:
    const double E_abs = 1e-8; 
    const double E_rel = 1e-8; 
    const double rho_max = 3;
    const double rho = 0.8;
    const double rho_min = 0.5;
    Vector_Function& vecFunc;
    std::vector<double> initial_value;
    double k;
    int p;
    int iteration_num;
    std::vector<std::vector<double>> results;
    std::string if_Richardson;
    std::vector<std::vector<double>> A;
    std::vector<double> b;
    std::vector<double> b_hat;
    std::vector<double> c;
};

class TimeIntegratorFactory 
{
public:
    static TimeIntegratorFactory& getInstance() 
    {
        static TimeIntegratorFactory instance;
        return instance;
    }

    template <typename... Args>
    std::unique_ptr<TimeIntegrator> createTimeIntegrator(const std::string& integratorName, Vector_Function& vecFunc, Args&&... args) {
        auto it = creators.find(integratorName);
        if (it != creators.end()) 
        {
            return it->second(vecFunc, std::make_tuple(std::forward<Args>(args)...));
        }
        return nullptr;
    }

private:
    TimeIntegratorFactory() 
    {
        creators["Adams_Bashforth"] = [](Vector_Function& vecFunc, const std::tuple<const std::vector<double>&, double, int, int, std::string>& args) 
        {
            auto& [initial_value, step_size, order, num_iterations, if_Richardson] = args;
            return std::make_unique<Adams_Bashforth>(vecFunc, initial_value, step_size, order, num_iterations, if_Richardson);
        };
        creators["Adams_Moulton"] = [](Vector_Function& vecFunc, const std::tuple<const std::vector<double>&, double, int, int, std::string>& args) 
        {
            auto& [initial_value, step_size, order, num_iterations, if_Richardson] = args;
            return std::make_unique<Adams_Moulton>(vecFunc, initial_value, step_size, order, num_iterations, if_Richardson);
        };
        creators["Backward_differentiation"] = [](Vector_Function& vecFunc, const std::tuple<const std::vector<double>&, double, int, int, std::string>& args) 
        {
            auto& [initial_value, step_size, order, num_iterations, if_Richardson] = args;
            return std::make_unique<Backward_differentiation>(vecFunc, initial_value, step_size, order, num_iterations, if_Richardson);
        };
        creators["Classical_RK"] = [](Vector_Function& vecFunc, const std::tuple<const std::vector<double>&, double, int, int, std::string>& args) 
        {
            auto& [initial_value, step_size, order, num_iterations, if_Richardson] = args;
            return std::make_unique<Classical_RK>(vecFunc, initial_value, step_size, order, num_iterations, if_Richardson);
        };
        creators["explicit_SDIRK"] = [](Vector_Function& vecFunc, const std::tuple<const std::vector<double>&, double, int, int, std::string>& args) 
        {
            auto& [initial_value, step_size, order, num_iterations, if_Richardson] = args;
            return std::make_unique<explicit_SDIRK>(vecFunc, initial_value, step_size, order, num_iterations, if_Richardson);
        };
        creators["Gauss_Legendre_RK"] = [](Vector_Function& vecFunc, const std::tuple<const std::vector<double>&, double, int, int, std::string>& args) 
        {
            auto& [initial_value, step_size, order, num_iterations, if_Richardson] = args;
            return std::make_unique<Gauss_Legendre_RK>(vecFunc, initial_value, step_size, order, num_iterations, if_Richardson);
        };
        creators["Fehlberg_embedded_RK"] = [](Vector_Function& vecFunc, const std::tuple<const std::vector<double>&, double, int, int, std::string>& args) 
        {
            auto& [initial_value, step_size, order, num_iterations, if_Richardson] = args;
            return std::make_unique<Fehlberg_embedded_RK>(vecFunc, initial_value, step_size, order, num_iterations, if_Richardson);
        };
        creators["Dormand_Prince_embedded_RK"] = [](Vector_Function& vecFunc, const std::tuple<const std::vector<double>&, double, int, int, std::string>& args) 
        {
            auto& [initial_value, step_size, order, num_iterations, if_Richardson] = args;
            return std::make_unique<Dormand_Prince_embedded_RK>(vecFunc, initial_value, step_size, order, num_iterations, if_Richardson);
        };
    }

    using Creator = std::function<std::unique_ptr<TimeIntegrator>(Vector_Function&, const std::tuple<const std::vector<double>&, double, int, int, std::string>&)>;
    std::map<std::string, Creator> creators;

    // Disallow copying and assignment
    TimeIntegratorFactory(const TimeIntegratorFactory&) = delete;
    TimeIntegratorFactory& operator=(const TimeIntegratorFactory&) = delete;
};

