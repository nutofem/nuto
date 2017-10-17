#include "BoostUnitTest.h"

#include <eigen3/Eigen/Core>
#include "mechanics/timeIntegration/RK4.h"

class HarmonicOscillator
{

    double m_gam;

public:
    HarmonicOscillator(double gam)
        : m_gam(gam)
    {
    }

    //! Equation for damped harmonic oscillator represented as
    //! first order system
    void operator()(const Eigen::Vector2d& w, Eigen::Vector2d& dxdt, double t)
    {
        dxdt[0] = w[1];
        dxdt[1] = -w[0] - 2. * m_gam * w[1];
    }

    //! Solution for damped harmonic oscillator at time t
    //! Initial conditions [displacement, velocity] w0
    double ExactSolution(const Eigen::Vector2d& w0, double t)
    {
        double om1 = sqrt(1. - m_gam * m_gam);
        double C1 = w0[0];
        double C2 = 1 / om1 * (w0[1] + m_gam * w0[0]);
        double result = exp(-m_gam * t) * (C1 * cos(om1 * t) + C2 * sin(om1 * t));
        return result;
    }
};


BOOST_AUTO_TEST_CASE(DampedHarmonicOscillator)
{
    NuTo::TimeIntegration::RK4<Eigen::Vector2d> ti;

    double stepSize = 0.01;
    int numSteps = 2000;

    HarmonicOscillator eq(0.15);
    Eigen::Vector2d y0 = {1.0, 0.0};
    Eigen::Vector2d y = y0;

    for (int i = 0; i < numSteps; i++)
    {
        double computed = y(0);
        double expected = eq.ExactSolution(y0, i * stepSize);
        BOOST_CHECK_CLOSE(computed, expected, 1e-3);
        y = ti.DoStep(eq, y, 0., stepSize);
    }
}
