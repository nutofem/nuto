#include "BoostUnitTest.h"

#include "math/shapes/Triangle.h"
#include "mechanics/integrationtypes/IntegrationTypeBase.h"
#include "mechanics/integrationtypes/IntegrationCompanion.h"

#include <vector>
#include <iostream>

//!@brief Integrates function f by quadrature
//!
//! The integration domain is given by the integration type.
//! @param f function to integrate
//! @param intType integrationtype to provide weights and integration points
//! @result computed integral: a scalar value
double integrate(std::function<double(Eigen::VectorXd)> f, NuTo::IntegrationTypeBase& intType)
{
    double result = 0.;
    for (int i = 0; i < intType.GetNumIntegrationPoints(); i++)
    {
        double y = f(intType.GetLocalIntegrationPointCoordinates(i));
        double w = intType.GetIntegrationPointWeight(i);
        result += w * y;
    }
    return (result);
}

double ExactIntegralMonomialTriangle(int order, int indx)
{
    std::vector<std::vector<double>> expectedResult = {
            {1. / 2},
            {1. / 6, 1. / 6},
            {1. / 12, 1. / 24, 1. / 12},
            {1. / 20, 1. / 60, 1. / 60, 1. / 20},
            {1. / 30, 1. / 120, 1. / 180, 1. / 120, 1. / 30},
            {1. / 42, 1. / 210, 1. / 420, 1. / 420, 1. / 210, 1. / 42},
            {1. / 56, 1. / 336, 1. / 840, 1. / 1120, 1. / 840, 1. / 336, 1. / 56},
            {1. / 72, 1. / 504, 1. / 1512, 1. / 2520, 1. / 2520, 1. / 1512, 1. / 504, 1. / 72},
            {1. / 90, 1. / 720, 1. / 2520, 1. / 5040, 1. / 6300, 1. / 5040, 1. / 2520, 1. / 720, 1. / 90},
            {1. / 110, 1. / 990, 1. / 3960, 1. / 9240, 1. / 13860, 1. / 13860, 1. / 9240, 1. / 3960, 1. / 990,
             1. / 110},
            {1. / 132, 1. / 1320, 1. / 5940, 1. / 15840, 1. / 27720, 1. / 33264, 1. / 27720, 1. / 15840, 1. / 5940,
             1. / 1320, 1. / 132},
            {1. / 156, 1. / 1716, 1. / 8580, 1. / 25740, 1. / 51480, 1. / 72072, 1. / 72072, 1. / 51480, 1. / 25740,
             1. / 8580, 1. / 1716, 1. / 156},
            {1. / 182, 1. / 2184, 1. / 12012, 1. / 40040, 1. / 90090, 1. / 144144, 1. / 168168, 1. / 144144, 1. / 90090,
             1. / 40040, 1. / 12012, 1. / 2184, 1. / 182},
            {1. / 210, 1. / 2730, 1. / 16380, 1. / 60060, 1. / 150150, 1. / 270270, 1. / 360360, 1. / 360360,
             1. / 270270, 1. / 150150, 1. / 60060, 1. / 16380, 1. / 2730, 1. / 210},
            {1. / 240, 1. / 3360, 1. / 21840, 1. / 87360, 1. / 240240, 1. / 480480, 1. / 720720, 1. / 823680,
             1. / 720720, 1. / 480480, 1. / 240240, 1. / 87360, 1. / 21840, 1. / 3360, 1. / 240},
            {1. / 272, 1. / 4080, 1. / 28560, 1. / 123760, 1. / 371280, 1. / 816816, 1. / 1361360, 1. / 1750320,
             1. / 1750320, 1. / 1361360, 1. / 816816, 1. / 371280, 1. / 123760, 1. / 28560, 1. / 4080, 1. / 272},
            {1. / 306, 1. / 4896, 1. / 36720, 1. / 171360, 1. / 556920, 1. / 1336608, 1. / 2450448, 1. / 3500640,
             1. / 3938220, 1. / 3500640, 1. / 2450448, 1. / 1336608, 1. / 556920, 1. / 171360, 1. / 36720, 1. / 4896,
             1. / 306},
            {1. / 342, 1. / 5814, 1. / 46512, 1. / 232560, 1. / 813960, 1. / 2116296, 1. / 4232592, 1. / 6651216,
             1. / 8314020, 1. / 8314020, 1. / 6651216, 1. / 4232592, 1. / 2116296, 1. / 813960, 1. / 232560, 1. / 46512,
             1. / 5814, 1. / 342},
            {1. / 380, 1. / 6840, 1. / 58140, 1. / 310080, 1. / 1162800, 1. / 3255840, 1. / 7054320, 1. / 12093120,
             1. / 16628040, 1. / 18475600, 1. / 16628040, 1. / 12093120, 1. / 7054320, 1. / 3255840, 1. / 1162800,
             1. / 310080, 1. / 58140, 1. / 6840, 1. / 380},
            {1. / 420,      1. / 7980,     1. / 71820,    1. / 406980,   1. / 1627920,  1. / 4883760,  1. / 11395440,
             1. / 21162960, 1. / 31744440, 1. / 38798760, 1. / 38798760, 1. / 31744440, 1. / 21162960, 1. / 11395440,
             1. / 4883760,  1. / 1627920,  1. / 406980,   1. / 71820,    1. / 7980,     1. / 420},
            {1. / 462,      1. / 9240,     1. / 87780,    1. / 526680,   1. / 2238390,  1. / 7162848,  1. / 17907120,
             1. / 35814240, 1. / 58198140, 1. / 77597520, 1. / 85357272, 1. / 77597520, 1. / 58198140, 1. / 35814240,
             1. / 17907120, 1. / 7162848,  1. / 2238390,  1. / 526680,   1. / 87780,    1. / 9240,     1. / 462},
            {1. / 506,       1. / 10626,     1. / 106260,    1. / 672980,    1. / 3028410,   1. / 10296594,
             1. / 27457584,  1. / 58837680,  1. / 102965940, 1. / 148728580, 1. / 178474296, 1. / 178474296,
             1. / 148728580, 1. / 102965940, 1. / 58837680,  1. / 27457584,  1. / 10296594,  1. / 3028410,
             1. / 672980,    1. / 106260,    1. / 10626,     1. / 506}};
    return expectedResult[order][indx];
}

void CheckTriangleIntegration(int polyOrder, int intTypeOrder)
{
    std::unique_ptr<NuTo::IntegrationTypeBase> intType = NuTo::CreateIntegrationType(NuTo::Triangle(), intTypeOrder);
    for (int n = 0; n <= polyOrder; n++)
    {
        int count = 0;
        for (int i = 0; i < n + 1; i++)
        {
            auto f = [n, i](Eigen::VectorXd x) { return (std::pow(x[0], i) * std::pow(x[1], n - i)); };
            double computedResult = integrate(f, *intType);
            BOOST_CHECK_CLOSE(computedResult, ExactIntegralMonomialTriangle(n, count), 1.e-11);
            count++;
        }
    }
}

BOOST_AUTO_TEST_CASE(IntegrateTriangle)
{
    int maxOrder = 19;
    for (int i = 1; i <= maxOrder; i++)
    {
        CheckTriangleIntegration(i, i);
        std::cout << "Done: order " << i << std::endl;
    }
}
