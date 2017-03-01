#include "BoostUnitTest.h"
#include <fakeit.hpp>
#include "mechanics/cell/CellData.h"

namespace NuTo
{
class Interpolation;
}


BOOST_AUTO_TEST_CASE(CellDataNodeValues)
{
    NuTo::NodeSimple n0 = NuTo::NodeSimple(Eigen::Vector2d({1, 1}));
    NuTo::NodeSimple n1 = NuTo::NodeSimple(Eigen::Vector2d({5, 1}));
    NuTo::NodeSimple n2 = NuTo::NodeSimple(Eigen::Vector2d({1, 7}));

    fakeit::Mock<NuTo::Interpolation> interpolation;

    NuTo::Element({&n0, &n1, &n2}, interpolation.get());
}
