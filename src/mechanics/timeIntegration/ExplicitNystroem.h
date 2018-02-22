#pragma once

#include <vector>

namespace NuTo
{
namespace TimeIntegration
{
template <typename TState>
class ExplicitNystroem
{

public:
    //! @brief Initialization with method specific parameters (butcher tableau)
    ExplicitNystroem(std::vector<std::vector<double>> aa, std::vector<double> bb1, std::vector<double> bb2,
                     std::vector<double> cc)
        : a(aa)
        , b1(bb1)
        , b2(bb2)
        , c(cc)
    {
    }


    //! @brief Performs one Nystroem step
    //! @param f A functor that returns the right hand side of the
    //! special second order differential equation w''=f(w,t)
    //!
    //! The signature of its call operator must be:
    //! operator()(const TState& w, TState& d2wdt2, double t)
    //! The return value is stored in d2wdt2
    //!
    //! @param w0 initial value
    //! @param v0 initial velocity
    //! @param t0 start time
    //! @param h step size (t-t0)
    //! @return pair of value and velocity after one Nystroem step
    template <typename TFunctor>
    std::pair<TState, TState> DoStep(TFunctor f, TState w0, TState v0, double t0, double h)
    {
        std::vector<TState> k(c.size(), w0);
        TState resultW = w0 + h * v0;
        TState resultV = v0;

        for (std::size_t i = 0; i < c.size(); i++)
        {
            double t = t0 + c[i] * h;
            TState wni = w0 + h * v0 * c[i];
            for (std::size_t j = 0; j < i; j++)
            {
                if (a[i][j] != 0)
                    wni += h * h * a[i][j] * k[j];
            }
            f(wni, k[i], t);
            resultW += h * h * b2[i] * k[i];
            resultV += h * b1[i] * k[i];
        }
        return std::make_pair(resultW, resultV);
    }


protected:
    std::vector<std::vector<double>> a;
    std::vector<double> b1;
    std::vector<double> b2;
    std::vector<double> c;
};
}
}
