#include "BoostUnitTest.h"
#include "TypeTraits.h"
#include "mechanics/interpolation/InterpolationTriangle.h"


/**
 * Define a InterpolationTest.h that runs these checks
 * for arbitrary interpolations
 */

BOOST_AUTO_TEST_CASE(InterpolationTriangleCopyMove)
{
    NuTo::Test::Copy<NuTo::InterpolationTriangle>();
    NuTo::Test::Move<NuTo::InterpolationTriangle>();
}

BOOST_AUTO_TEST_CASE(InterpolationTriangleN)
{
    NuTo::InterpolationTriangle interpolation(NuTo::eInterpolation::GAUSS, 1);

    for (int iNode = 0; iNode < interpolation.GetNumNodes(); ++iNode)
    {
        auto localNodeCoordinates = interpolation.GetLocalCoords(iNode);
        auto N                    = interpolation.GetShapeFunctions(localNodeCoordinates);
        for (int i = 0; i < interpolation.GetNumNodes(); ++i)
        {
            if (i == iNode)
                BOOST_CHECK_CLOSE(N[i], 1, 1.e-10);
            else
                BOOST_CHECK_SMALL(N[i], 1.e-10);
        }
    }
}

BOOST_AUTO_TEST_CASE(InterpolationTriangleB)
{
    // check derviatives via CDF.
}
