#include "BoostUnitTest.h"
#include <type_traits>

#include "mechanics/elements/Element.h"
#include "mechanics/interpolation/InterpolationTriangle.h"

BOOST_AUTO_TEST_CASE(ElementCopyMove)
{
    BOOST_CHECK(std::is_copy_constructible<NuTo::Element>::value);
    BOOST_CHECK(std::is_move_constructible<NuTo::Element>::value);
}


struct TestElement : public NuTo::Element
{
    TestElement()
        : NuTo::Element({&n0, &n1, &n2}, interpolation)
    {
    }

    NuTo::NodeSimple n0                       = NuTo::NodeSimple(Eigen::Vector2d({1, 1}));
    NuTo::NodeSimple n1                       = NuTo::NodeSimple(Eigen::Vector2d({5, 1}));
    NuTo::NodeSimple n2                       = NuTo::NodeSimple(Eigen::Vector2d({1, 7}));
    NuTo::InterpolationTriangle interpolation = NuTo::InterpolationTriangle(NuTo::eInterpolation::GAUSS, 1, 2);
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
