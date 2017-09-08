#include "BoostUnitTest.h"
#include <boost/test/output_test_stream.hpp>
#include <type_traits>

#include "mechanics/interpolation/CellInterpolationIGA.h"
#include <iostream>

BOOST_AUTO_TEST_CASE(ElementCopyMove)
{
    BOOST_CHECK(std::is_copy_constructible<NuTo::CellInterpolationIGA<1>>::value);
    BOOST_CHECK(std::is_move_constructible<NuTo::CellInterpolationIGA<1>>::value);

    BOOST_CHECK(std::is_copy_constructible<NuTo::CellInterpolationIGA<2>>::value);
    BOOST_CHECK(std::is_move_constructible<NuTo::CellInterpolationIGA<2>>::value);
}


BOOST_AUTO_TEST_CASE(ExtractNodeValues)
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

    std::array<int, 1> knotIDsCell = {4};
    NuTo::CellInterpolationIGA<1> iga(knotIDsCell, curve);
    NuTo::NodeValues nodeValues = iga.ExtractNodeValues();

    // ip coordinates are passed in
    Eigen::VectorXd param(1); param << 0.;
    Eigen::VectorXd shapefuns = iga.GetShapeFunctions(param);

    Eigen::Vector2d result(0,0);

    result(0) = shapefuns(0)*nodeValues(0) + shapefuns(1)*nodeValues(2) + shapefuns(2)*nodeValues(4);
    result(1) = shapefuns(0)*nodeValues(1) + shapefuns(1)*nodeValues(3) + shapefuns(2)*nodeValues(5);

    // curve parameter coordinates are passed in
    param << 0.375;
    BoostUnitTest::CheckVector(result, curve.Evaluate(param), 2);
}

