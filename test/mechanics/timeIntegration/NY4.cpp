#include "BoostUnitTest.h"

#include <iostream>
#include <Eigen/Core>
#include "mechanics/timeIntegration/NY4.h"

class HarmonicOscillatorUndamped
{

public:
    HarmonicOscillatorUndamped()
        : m_gam(gam)
    {
    }

    //! Equation for harmonic oscillator represented as
    //! second order system
    void operator()(const double& w, double& d2wdt2, double t)
    {
        d2wdt2 = -w;
    }

    //! Solution for undamped harmonic oscillator at time t
    //! Initial conditions [displacement, velocity] w0,v0
    double exactSolution(const double& w0, const double& v0, double t)
    {
        return w0 * cos(t) + v0 * sin(t);
    }
};


BOOST_AUTO_TEST_CASE(HarmonicOscillator)
{
    NuTo::TimeIntegration::NY4<double> ti;

    double stepSize = 0.01;
    int numSteps = 2000;

    HarmonicOscillatorSecondOrder eq(0.15);
    double y0 = 1.0;
    double v0 = 0.0;
    std::pair<double, double> y = std::make_pair(y0, v0);

    for (int i = 0; i < numSteps; i++)
    {
        double computed = y.first;
        double expected = eq.exactSolution(y0, v0, i * stepSize);
        BOOST_CHECK_SMALL(computed - expected, 1e-4);
        y = ti.DoStep(eq, y.first, y.second, i * stepSize, stepSize);
    }
}
