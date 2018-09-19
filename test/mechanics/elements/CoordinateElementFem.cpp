#include "BoostUnitTest.h"
#include <type_traits>

#include "nuto/mechanics/elements/CoordinateElementFem.h"
#include "nuto/mechanics/interpolation/InterpolationTriangleLinear.h"


BOOST_AUTO_TEST_CASE(ElementCopyMove)
{
    BOOST_CHECK(std::is_copy_constructible<NuTo::CoordinateElementFem>::value);
    BOOST_CHECK(std::is_move_constructible<NuTo::CoordinateElementFem>::value);
}


NuTo::CoordinateNode n0 = NuTo::CoordinateNode(Eigen::Vector2d({1, 1}));
NuTo::CoordinateNode n1 = NuTo::CoordinateNode(Eigen::Vector2d({5, 1}));
NuTo::CoordinateNode n2 = NuTo::CoordinateNode(Eigen::Vector2d({1, 7}));
NuTo::InterpolationTriangleLinear interpolation;

NuTo::CoordinateElementFem TestElement()
{
    return NuTo::CoordinateElementFem({n0, n1, n2}, interpolation);
}

BOOST_AUTO_TEST_CASE(ExtractNodeValues)
{
    Eigen::VectorXd nodeValues = TestElement().ExtractCoordinates();
    BoostUnitTest::CheckVector(nodeValues, std::vector<double>{1, 1, 5, 1, 1, 7}, 6);
}

BOOST_AUTO_TEST_CASE(Interpolation)
{
    Eigen::VectorXd nodeValues = TestElement().ExtractCoordinates();
    Eigen::MatrixXd N = TestElement().GetNMatrix(Eigen::Vector2d(0.5, 0.5));
    Eigen::Vector2d interpolatedValues = N * nodeValues;
    BoostUnitTest::CheckVector(interpolatedValues, std::vector<double>{3, 4}, 2);
}
