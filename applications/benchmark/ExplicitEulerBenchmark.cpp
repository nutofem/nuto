#include "Benchmark.h"
#include <cmath>
#include <boost/numeric/odeint/stepper/runge_kutta4_classic.hpp>
#include <array>
#include "base/MakeArray.h"

/**
 *
 * Solving y' = y with y(0) = 1
 * exact solution: exp(t)
 *
 */

constexpr double dt = 0.001;
constexpr double tEnd = 1;
constexpr double tolerance = 1.e-2;

void f(double y, double& dydt, const double /* t */)
{
    dydt = y;
}

void Check(double y)
{
    if (std::abs(std::exp(tEnd) - y) > tolerance)
        throw;
}

BENCHMARK(Euler, odeint, runner)
{
    while (runner.KeepRunningIterations(1000))
    {
        boost::numeric::odeint::runge_kutta4_classic<double> euler;
        double y = 1;
        for (double t = 0; t < tEnd; t += dt)
            euler.do_step(f, y, t, dt);
        Check(y);
    }
}
//
// BENCHMARK(Euler, hardcode, runner)
//{
//    while (runner.KeepRunningIterations(1000))
//    {
//        double y = 1;
//        for (double t = 0; t < tEnd; t += dt)
//            y += dt * y;
//        Check(y);
//    }
//}

template <typename TState, typename TButcher>
class ExplicitRungeKutta
{

public:
    //! @brief Performs one RungeKutta step
    //! @param f A functor that returns the right hand side of the differential equation
    //!
    //! The signature of its call operator must be:
    //! operator()(const TState& w, TState& dwdt, double t)
    //! The return value is stored in dwdt
    //!
    //! @param w0 initial value
    //! @param t0 start time
    //! @param h step size (t-t0)
    //! @return value after one RungeKutta step
    template <typename F>
    TState DoStep(F f, TState w0, double t0, double h)
    {
        std::array<TState, TButcher::Num> k = NuTo::make_array_n<TButcher::Num>(w0);
        TState result = w0;

        for (std::size_t i = 0; i < TButcher::Num; i++)
        {
            double t = t0 + TButcher::C[i] * h;
            TState wni = w0;
            for (std::size_t j = 0; j < i; j++)
            {
                if (TButcher::A[i][j] != 0)
                    wni += h * TButcher::A[i][j] * k[j];
            }
            f(wni, k[i], t);
            result += h * TButcher::B[i] * k[i];
        }
        return result;
    }
};


struct BtRungeKutta4
{
    static constexpr std::array<std::array<double, 4>, 4> A = {{{0, 0, 0, 0}, {0.5, 0, 0, 0}, {0.0, 0.5, 0, 0}, {0., 0., 1., 0}}};
    static constexpr std::array<double, 4> B = {1. / 6., 1. / 3., 1. / 3., 1. / 6.};
    static constexpr std::array<double, 4> C = {1. / 6., 1. / 3., 1. / 3., 1. / 6.};
    static constexpr int Num = 4;
};

constexpr std::array<std::array<double, 4>, 4> BtRungeKutta4::A;
constexpr std::array<double, 4> BtRungeKutta4::B;
constexpr std::array<double, 4> BtRungeKutta4::C;

using ExplicitEuler = ExplicitRungeKutta<double, BtRungeKutta4>;

BENCHMARK(Euler, NuToRK, runner)
{
    while (runner.KeepRunningIterations(1000))
    {
        ExplicitEuler euler;
        double y = 1;
        for (double t = 0; t < tEnd; t += dt)
            y = euler.DoStep(f, y, t, dt);
        Check(y);
    }
}
