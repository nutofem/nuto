#include "BoostUnitTest.h"
#include "InterpolationTests.h"
#include "mechanics/interpolation/InterpolationBrickLobatto.h"

#include <iostream>

std::vector<Eigen::VectorXd> GetTestPoints()
{
    constexpr double a = 0.577350269189626;
    return {Eigen::Vector3d({-a, -a, -a}), Eigen::Vector3d({-a, a, -a}), Eigen::Vector3d({a, a, -a}),
            Eigen::Vector3d({a, -a, -a}),  Eigen::Vector3d({-a, -a, a}), Eigen::Vector3d({-a, a, a}),
            Eigen::Vector3d({a, a, a}),    Eigen::Vector3d({a, -a, a})};
}

BOOST_AUTO_TEST_CASE(InterpolationBrickLobatto)
{
    int order = 4;
    NuTo::InterpolationBrickLobatto interpolation(1, order);

    std::vector<Eigen::VectorXd> rPoints = GetTestPoints();

    BOOST_TEST_MESSAGE("Checking copy move...");
    std::cout << "Checking copy move..." << std::endl;
    NuTo::Test::CheckCopyMove<NuTo::InterpolationBrickLobatto>();
    BOOST_TEST_MESSAGE("Checking match of shape functions and node positions...");
    std::cout << "Shape / Node" << std::endl;
    NuTo::Test::CheckShapeFunctionsAndNodePositions(interpolation);
    BOOST_TEST_MESSAGE("Checking partition of unity...");
    std::cout << "PoU" << std::endl;
    NuTo::Test::CheckPartitionOfUnity(interpolation, rPoints);
    BOOST_TEST_MESSAGE("Checking shape function derivatives via CDF...");
    std::cout << "Derivatives" << std::endl;
    NuTo::Test::CheckDerivativeShapeFunctionsCDF(interpolation, rPoints);
}
