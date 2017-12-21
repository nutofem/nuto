#include "BoostUnitTest.h"
#include <type_traits>

#include "mechanics/elements/ElementFem.h"
#include "mechanics/interpolation/InterpolationTriangleLinear.h"


BOOST_AUTO_TEST_CASE(ElementCopyMove)
{
    BOOST_CHECK(std::is_copy_constructible<NuTo::ElementFem>::value);
    BOOST_CHECK(std::is_move_constructible<NuTo::ElementFem>::value);
}


NuTo::NodeSimple n0 = NuTo::NodeSimple(Eigen::Vector2d({1, 1}));
NuTo::NodeSimple n1 = NuTo::NodeSimple(Eigen::Vector2d({5, 1}));
NuTo::NodeSimple n2 = NuTo::NodeSimple(Eigen::Vector2d({1, 7}));
NuTo::InterpolationTriangleLinear interpolation;

NuTo::ElementFem TestElement()
{
    return NuTo::ElementFem({n0, n1, n2}, interpolation);
}

BOOST_AUTO_TEST_CASE(ExtractNodeValues)
{
    NuTo::NodeValues nodeValues = TestElement().ExtractNodeValues();
    BoostUnitTest::CheckVector(nodeValues, std::vector<double>{1, 1, 5, 1, 1, 7}, 6);
}

BOOST_AUTO_TEST_CASE(Interpolation)
{
    NuTo::NodeValues nodeValues = TestElement().ExtractNodeValues();
    NuTo::NMatrix N = TestElement().GetNMatrix(Eigen::Vector2d(0.5, 0.5));
    Eigen::Vector2d interpolatedValues = N * nodeValues;
    BoostUnitTest::CheckVector(interpolatedValues, std::vector<double>{3, 4}, 2);
}
