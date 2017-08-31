#include "BoostUnitTest.h"
#include <type_traits>

#include "mechanics/IGA/NURBS.h"

BOOST_AUTO_TEST_CASE(NURBSCopyMove)
{
    BOOST_CHECK(std::is_copy_constructible<NuTo::NURBS<1>>::value);
    BOOST_CHECK(std::is_move_constructible<NuTo::NURBS<1>>::value);

    BOOST_CHECK(std::is_copy_constructible<NuTo::NURBS<2>>::value);
    BOOST_CHECK(std::is_move_constructible<NuTo::NURBS<2>>::value);
}

BOOST_AUTO_TEST_CASE(NURBS)
{
    // IGA geometry of a circle (NURBS curve should exactly fit the circle)
    std::vector<std::vector<NuTo::NodeSimple>> controlPoints;
    std::vector<NuTo::NodeSimple> row1;

    row1.push_back(NuTo::NodeSimple(Eigen::Vector2d({0,  -1})));
    row1.push_back(NuTo::NodeSimple(Eigen::Vector2d({-1, -1})));
    row1.push_back(NuTo::NodeSimple(Eigen::Vector2d({-1,  0})));
    row1.push_back(NuTo::NodeSimple(Eigen::Vector2d({-1,  1})));
    row1.push_back(NuTo::NodeSimple(Eigen::Vector2d({0,   1})));
    row1.push_back(NuTo::NodeSimple(Eigen::Vector2d({1,   1})));
    row1.push_back(NuTo::NodeSimple(Eigen::Vector2d({1,   0})));
    row1.push_back(NuTo::NodeSimple(Eigen::Vector2d({1,  -1})));
    row1.push_back(NuTo::NodeSimple(Eigen::Vector2d({0,  -1})));

    controlPoints.push_back(row1);

    std::vector<double> knots1D = {0, 0, 0, 1/4., 1/4., 1/2., 1/2., 3/4., 3/4., 1, 1, 1};
    std::array<std::vector<double>, 1> knots = {knots1D};

    std::vector<double> weights1D = {1, std::sqrt(2)/2, 1, std::sqrt(2)/2, 1, std::sqrt(2)/2, 1, std::sqrt(2)/2, 1};
    std::vector<std::vector<double>> weights = {weights1D};

    std::array<int, 1> degree = {2};

    NuTo::NURBS<1> curve(knots, controlPoints, weights, degree);

    Eigen::Matrix<double, 1, 1> vec;

    for(int i = 0; i < 100; i++)
    {
        vec << i/100.;
        Eigen::VectorXd vel = curve.Evaluate(vec);
        BOOST_CHECK_CLOSE(vel(0)*vel(0)+vel(1)*vel(1), 1 , 1.e-10);
    }
}
