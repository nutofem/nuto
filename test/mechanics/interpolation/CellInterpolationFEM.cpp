#include "BoostUnitTest.h"
#include <type_traits>

#include "mechanics/interpolation/CellInterpolationFEM.h"
#include "mechanics/interpolation/InterpolationTriangleLinear.h"


BOOST_AUTO_TEST_CASE(ElementCopyMove)
{
    BOOST_CHECK(std::is_copy_constructible<NuTo::CellInterpolationFEM>::value);
    BOOST_CHECK(std::is_move_constructible<NuTo::CellInterpolationFEM>::value);
}


struct TestElement : public NuTo::CellInterpolationFEM
{
    TestElement()
        : NuTo::CellInterpolationFEM({&n0, &n1, &n2}, interpolation)
    {
    }

    NuTo::NodeSimple n0 = NuTo::NodeSimple(Eigen::Vector2d({1, 1}));
    NuTo::NodeSimple n1 = NuTo::NodeSimple(Eigen::Vector2d({5, 1}));
    NuTo::NodeSimple n2 = NuTo::NodeSimple(Eigen::Vector2d({1, 7}));
    NuTo::InterpolationTriangleLinear interpolation = NuTo::InterpolationTriangleLinear(2);
};

BOOST_AUTO_TEST_CASE(ExtractNodeValues)
{
    NuTo::NodeValues nodeValues = TestElement().ExtractNodeValues();
    BoostUnitTest::CheckVector(nodeValues, std::vector<double>{1, 1, 5, 1, 1, 7}, 6);
}

