#include "BoostUnitTest.h"

#include <iostream>
#include <Eigen/Core>
#include "mechanics/timeIntegration/NY4.h"
#include "mechanics/timeIntegration/NY5.h"
#include "mechanics/timeIntegration/NYVelocityVerlet.h"

class HarmonicOscillatorUndamped
{

public:
    HarmonicOscillatorUndamped()
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

void Run(NuTo::TimeIntegration::ExplicitNystroem<double> ti)
{
    double stepSize = 0.01;
    int numSteps = 2000;

    HarmonicOscillatorUndamped eq;
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

BOOST_AUTO_TEST_CASE(NY4)
{
    Run(NuTo::TimeIntegration::NY4<double>());
}

BOOST_AUTO_TEST_CASE(NY5)
{
    Run(NuTo::TimeIntegration::NY5<double>());
}

BOOST_AUTO_TEST_CASE(NYVelocityVerlet)
{
    Run(NuTo::TimeIntegration::NYVelocityVerlet<double>());
}
