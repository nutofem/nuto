#include "BoostUnitTest.h"
#include <type_traits>

#include "mechanics/elements/ElementSimple.h"
#include "mechanics/interpolation/InterpolationTriangleLinear.h"

BOOST_AUTO_TEST_CASE(ElementCopyMove)
{
    BOOST_CHECK(std::is_copy_constructible<NuTo::ElementSimple>::value);
    BOOST_CHECK(std::is_move_constructible<NuTo::ElementSimple>::value);
}


struct TestElement : public NuTo::ElementSimple
{
    TestElement()
        : NuTo::ElementSimple({&n0, &n1, &n2}, interpolation)
    {
    }

    NuTo::NodeSimple n0                       = NuTo::NodeSimple(Eigen::Vector2d({1, 1}));
    NuTo::NodeSimple n1                       = NuTo::NodeSimple(Eigen::Vector2d({5, 1}));
    NuTo::NodeSimple n2                       = NuTo::NodeSimple(Eigen::Vector2d({1, 7}));
    NuTo::InterpolationTriangleLinear interpolation = NuTo::InterpolationTriangleLinear(2);
};

BOOST_AUTO_TEST_CASE(ElementExtractNodeValues)
{
    NuTo::NodeValues nodeValues = TestElement().ExtractNodeValues();
    BoostUnitTest::CheckVector(nodeValues, std::vector<double>{1, 1, 5, 1, 1, 7}, 6);
}


BOOST_AUTO_TEST_CASE(ElementInterpolate)
{
    auto interpolatedValues     = TestElement().Interpolate(Eigen::Vector2d({0.5, 0.5}));
    std::vector<double> correct = {3., 4.};
    BoostUnitTest::CheckVector(interpolatedValues, correct, 2);
}
