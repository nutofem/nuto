#include "Eigen/Core"
#include <boost/numeric/odeint.hpp>
#include "nuto/math/EigenOdeintCompatibility.h"

/* Simple first order system (two harmonic oscillators)
 *
 * u1 = cos omega t, u2 = sin omega t
 * u1' = -omega u2, u2' = u1
 *
*/
class Vibrations
{
public:
    Vibrations(double om)
        : mOmega(om)
    {
    }

    void operator()(const Eigen::VectorXd& x, Eigen::VectorXd& dxdt, const double)
    {
        dxdt(0) = -mOmega * x(1);
        dxdt(1) = mOmega * x(0);
    }

private:
    double mOmega;
};

class Observer
{
public:
    Observer(std::vector<Eigen::VectorXd>& states, std::vector<double>& times)
        : m_states(states)
        , m_times(times)
    {
    }

    void operator()(const Eigen::VectorXd& x, double t)
    {
        m_states.push_back(x);
        m_times.push_back(t);
    }

private:
    std::vector<Eigen::VectorXd>& m_states;
    std::vector<double>& m_times;
};

int main()
{
    using namespace boost::numeric::odeint;

    Eigen::VectorXd x(2);
    x(0) = 1.;
    x(1) = 0.;

    Vibrations eq(2 * M_PI);

    std::vector<Eigen::VectorXd> x_vec;
    std::vector<double> times;

    runge_kutta_cash_karp54<Eigen::VectorXd> stepper;

    double t0 = 0.;
    double tF = 5.;
    double dt = 0.1; // initial step size

    integrate_adaptive(stepper, eq, x, t0, tF, dt, Observer(x_vec, times));

    std::cout << std::fixed;
    std::cout << "Time      u1          u2\n";
    for (size_t i = 0; i < times.size(); i++)
    {
        std::cout << std::left << std::setw(6) << std::setprecision(1) << times[i] << std::right << std::setw(12)
                  << std::setprecision(6) << x_vec[i][0] << std::setw(12) << x_vec[i][1] << '\n';
    }
}
