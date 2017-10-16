#pragma once

#include <vector>

namespace NuTo
{
namespace TimeIntegration
{
template <typename state>
class ExplicitRungeKutta
{

public:
    //! @brief Initialization with method specific parameters (butcher tableau)
    ExplicitRungeKutta(std::vector<std::vector<double>> aa, std::vector<double> bb, std::vector<double> cc)
    {
        for (auto elm : aa)
        {
            a.push_back(elm);
        }
        for (auto elm : bb)
        {
            b.push_back(elm);
        }
        for (auto elm : cc)
        {
            c.push_back(elm);
        }
    }

    //! @brief Performs one RungeKutta step
    template <typename F>
    state doStep(F f, state w0, double t0, double h)
    {
        std::vector<state> k;
        for (size_t i = 0; i < c.size(); i++)
        {
            k.push_back(state(w0));
        }
        state result = w0;
        state wni = w0;
        double t;

        for (std::size_t i = 0; i < c.size(); i++)
        {
            t = t0 + c[i] * h;
            wni = w0;
            for (std::size_t j = 0; j < i; j++)
            {
                wni += h * a[i][j] * k[j];
            }
            f(wni, k[i], t);
            result += h * b[i] * k[i];
        }
        return result;
    }

protected:
    std::vector<std::vector<double>> a;
    std::vector<double> b;
    std::vector<double> c;
};
}
}
