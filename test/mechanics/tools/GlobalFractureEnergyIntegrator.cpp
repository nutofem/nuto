#include "BoostUnitTest.h"
#include "nuto/mechanics/tools/GlobalFractureEnergyIntegrator.h"

BOOST_AUTO_TEST_CASE(IntegrateFunction)
{
    int num = 10000;
    Eigen::VectorXd displ(num);
    Eigen::VectorXd force(num);
    for (int i = 0; i < num; ++i)
    {
        double x = 2. / (num - 1) * i;
        displ[i] = x;
        force[i] = 1 - (x - 1) * (x - 1);
    }
    NuTo::Tools::GlobalFractureEnergyIntegrator integrator(force, displ);
    // integrating to zero
    BOOST_CHECK_CLOSE(integrator.IntegrateSofteningCurve(1, 0), 4. / 3., 1.e-5);
    BOOST_CHECK_CLOSE(integrator.IntegrateSofteningCurve(3, 0), 4. / 9., 1.e-5);
    BOOST_CHECK_CLOSE(integrator.IntegrateSofteningCurve(3), 4. / 9., 1.e-5);

    // integrating to threshold f = 0.75 after reaching the peak load. (Integrate [0 .. 1.5])
    BOOST_CHECK_CLOSE(integrator.IntegrateSofteningCurve(1, 0.75), 1.125, 1.e-5);
}

//! this may correspond to a load displacement curve in compression
BOOST_AUTO_TEST_CASE(IntegrateFunctionNegative)
{
    int num = 10000;
    Eigen::VectorXd displ(num);
    Eigen::VectorXd force(num);
    for (int i = 0; i < num; ++i)
    {
        double x = 2. / (num - 1) * i;
        displ[i] = -x;
        force[i] = -(1 - (x - 1) * (x - 1));
    }

    NuTo::Tools::GlobalFractureEnergyIntegrator integrator(force, displ);
    // integrating to zero
    BOOST_CHECK_CLOSE(integrator.IntegrateSofteningCurve(1, 0), 4. / 3., 1.e-5);
    BOOST_CHECK_CLOSE(integrator.IntegrateSofteningCurve(3, 0), 4. / 9., 1.e-5);
    BOOST_CHECK_CLOSE(integrator.IntegrateSofteningCurve(3), 4. / 9., 1.e-5);

    // integrating to threshold f = 0.75 after reaching the peak load. (Integrate [0 .. 1.5])
    BOOST_CHECK_CLOSE(integrator.IntegrateSofteningCurve(1, 0.75), 1.125, 1.e-5);
}
