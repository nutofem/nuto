#include "BoostUnitTest.h"
#include "mechanics/mesh/Mesh.h"
#include "mechanics/interpolation/InterpolationTriangleLinear.h"
#include <iostream>

BOOST_AUTO_TEST_CASE(MeshBasics)
{
    NuTo::Mesh mesh;
    NuTo::InterpolationTriangleLinear interpolationTriangle(2);
    auto& interpolation = mesh.CreateInterpolation(interpolationTriangle);

    BOOST_CHECK_EQUAL(interpolation.GetNumNodes(), 3);
    BOOST_CHECK_EQUAL(interpolation.GetDofDimension(), 2);


    auto& n0 = mesh.CreateNode(Eigen::Vector2d({1,0}));
    auto& n1 = mesh.CreateNode(Eigen::Vector2d({2,0}));
    auto& n2 = mesh.CreateNode(Eigen::Vector2d({0,3}));

    BOOST_CHECK_EQUAL(n0.GetNumValues(), 2);
    BOOST_CHECK_EQUAL(n1.GetNumValues(), 2);
    BOOST_CHECK_EQUAL(n2.GetNumValues(), 2);

    auto& e0 = mesh.CreateElement({&n0, &n1, &n2}, interpolation);

    BoostUnitTest::CheckVector(e0.ExtractNodeValues(), std::vector<double>({1,0,2,0,0,3}), 6);
}
