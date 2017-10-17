#pragma once

#include <vector>

namespace NuTo
{
namespace TimeIntegration
{
template <typename Tstate>
class ExplicitRungeKutta
{

public:
    //! @brief Initialization with method specific parameters (butcher tableau)
    ExplicitRungeKutta(std::vector<std::vector<double>> aa, std::vector<double> bb, std::vector<double> cc)
        : a(aa)
        , b(bb)
        , c(cc)
    {
    }

    //! @brief Performs one RungeKutta step
    //! @param f A functor that returns the right hand side of the differential equation
    //!
    //! The signature of its call operator must be:
    //! operator()(const Tstate& w, Tstate& dwdt, double t)
    //! The return value is stored in dwdt
    //!
    //! @param w0 initial value
    //! @param t0 start time
    //! @param h step size (t-t0)
    //! @return value after one RungeKutta step
    template <typename F>
    Tstate DoStep(F f, Tstate w0, double t0, double h)
    {
        std::vector<Tstate> k(c.size(), w0);
        Tstate result = w0;

        for (std::size_t i = 0; i < c.size(); i++)
        {
            double t = t0 + c[i] * h;
            Tstate wni = w0;
            for (std::size_t j = 0; j < i; j++)
            {
                if (a[i][j] != 0)
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
