#include "BoostUnitTest.h"
#include <fakeit.hpp>
#include "mechanics/cell/Cell.h"
#include "mechanics/elements/Element.h"
#include "mechanics/interpolation/InterpolationTriangle.h"

BOOST_AUTO_TEST_CASE(CellLetsSee)
{
    //fakeit::Mock<NuTo::Interpolation> mockInterpolation;
    NuTo::InterpolationTriangle interpolation(NuTo::eInterpolation::GAUSS, 1);

    NuTo::NodeSimple n0(Eigen::Vector2d({1, 1}));
    NuTo::NodeSimple n1(Eigen::Vector2d({2, 1}));
    NuTo::NodeSimple n2(Eigen::Vector2d({1, 4}));

    //NuTo::Element coordinateElement({&n0, &n1, &n2}, mockInterpolation.get());
    NuTo::Element coordinateElement({&n0, &n1, &n2}, interpolation);


    NuTo::Cell cell(coordinateElement);
}
