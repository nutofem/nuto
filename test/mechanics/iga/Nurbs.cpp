#include "BoostUnitTest.h"
#include <type_traits>
#include <iostream>

#include "nuto/mechanics/iga/Nurbs.h"

BOOST_AUTO_TEST_CASE(NurbsCopyMove)
{
    BOOST_CHECK(std::is_copy_constructible<NuTo::Nurbs<1>>::value);
    BOOST_CHECK(std::is_move_constructible<NuTo::Nurbs<1>>::value);

    BOOST_CHECK(std::is_copy_constructible<NuTo::Nurbs<2>>::value);
    BOOST_CHECK(std::is_move_constructible<NuTo::Nurbs<2>>::value);
}

BOOST_AUTO_TEST_CASE(NurbsCurve)
{
    // IGA geometry of a circle (Nurbs curve should exactly fit the circle)
    std::vector<std::vector<NuTo::DofNode*>> controlPoints;
    std::vector<NuTo::DofNode*> row1;

    NuTo::DofNode n0 = NuTo::DofNode(Eigen::Vector2d({0, -1}));
    row1.push_back(&n0);

    NuTo::DofNode n1 = NuTo::DofNode(Eigen::Vector2d({-1, -1}));
    row1.push_back(&n1);

    NuTo::DofNode n2 = NuTo::DofNode(Eigen::Vector2d({-1, 0}));
    row1.push_back(&n2);

    NuTo::DofNode n3 = NuTo::DofNode(Eigen::Vector2d({-1, 1}));
    row1.push_back(&n3);

    NuTo::DofNode n4 = NuTo::DofNode(Eigen::Vector2d({0, 1}));
    row1.push_back(&n4);

    NuTo::DofNode n5 = NuTo::DofNode(Eigen::Vector2d({1, 1}));
    row1.push_back(&n5);

    NuTo::DofNode n6 = NuTo::DofNode(Eigen::Vector2d({1, 0}));
    row1.push_back(&n6);

    NuTo::DofNode n7 = NuTo::DofNode(Eigen::Vector2d({1, -1}));
    row1.push_back(&n7);

    NuTo::DofNode n8 = NuTo::DofNode(Eigen::Vector2d({0, -1}));
    row1.push_back(&n8);

    controlPoints.push_back(row1);

    std::vector<double> knots1D = {0, 0, 0, 1 / 4., 1 / 4., 1 / 2., 1 / 2., 3 / 4., 3 / 4., 1, 1, 1};
    std::array<std::vector<double>, 1> knots = {knots1D};

    std::vector<double> weights1D = {1, std::sqrt(2) / 2, 1, std::sqrt(2) / 2, 1, std::sqrt(2) / 2,
                                     1, std::sqrt(2) / 2, 1};
    std::vector<std::vector<double>> weights = {weights1D};

    std::array<int, 1> degree = {2};

    NuTo::Nurbs<1> curve(knots, controlPoints, weights, degree);

    Eigen::Matrix<double, 1, 1> vec;

    for (int i = 0; i < 100; i++)
    {
        vec << i / 100.;
        Eigen::VectorXd vel = curve.Evaluate(vec);
        BOOST_CHECK_CLOSE(vel(0) * vel(0) + vel(1) * vel(1), 1, 1.e-10);
    }

    Eigen::Vector2d derivativeCheck(-1, 0);
    vec << 0;
    Eigen::VectorXd derivative = curve.Evaluate(vec, 1);
    derivative.normalize();
    BoostUnitTest::CheckEigenMatrix(derivative, derivativeCheck);

    vec << 0.25;
    derivativeCheck << 0, 1;
    derivative = curve.Evaluate(vec, 1);
    derivative.normalize();
    BoostUnitTest::CheckEigenMatrix(derivative, derivativeCheck);

    vec << 0.5;
    derivativeCheck << 1, 0;
    derivative = curve.Evaluate(vec, 1);
    derivative.normalize();
    BoostUnitTest::CheckEigenMatrix(derivative, derivativeCheck);

    vec << 0.75;
    derivativeCheck << 0, -1;
    derivative = curve.Evaluate(vec, 1);
    derivative.normalize();
    BoostUnitTest::CheckEigenMatrix(derivative, derivativeCheck);

    vec << 1;
    derivativeCheck << -1, 0;
    derivative = curve.Evaluate(vec, 1);
    derivative.normalize();
    BoostUnitTest::CheckEigenMatrix(derivative, derivativeCheck);
}

BOOST_AUTO_TEST_CASE(NurbsSurface)
{

    // IGA geometry of a circle (Nurbs curve should exactly fit the circle)
    std::vector<std::vector<NuTo::DofNode*>> controlPoints;

    std::vector<NuTo::DofNode*> row1;
    NuTo::DofNode n0 = NuTo::DofNode(Eigen::Vector3d({0, -1, 0}));
    row1.push_back(&n0);
    NuTo::DofNode n1 = NuTo::DofNode(Eigen::Vector3d({-1, -1, 0}));
    row1.push_back(&n1);
    NuTo::DofNode n2 = NuTo::DofNode(Eigen::Vector3d({-1, 0, 0}));
    row1.push_back(&n2);
    NuTo::DofNode n3 = NuTo::DofNode(Eigen::Vector3d({-1, 1, 0}));
    row1.push_back(&n3);
    NuTo::DofNode n4 = NuTo::DofNode(Eigen::Vector3d({0, 1, 0}));
    row1.push_back(&n4);
    NuTo::DofNode n5 = NuTo::DofNode(Eigen::Vector3d({1, 1, 0}));
    row1.push_back(&n5);
    NuTo::DofNode n6 = NuTo::DofNode(Eigen::Vector3d({1, 0, 0}));
    row1.push_back(&n6);
    NuTo::DofNode n7 = NuTo::DofNode(Eigen::Vector3d({1, -1, 0}));
    row1.push_back(&n7);
    NuTo::DofNode n8 = NuTo::DofNode(Eigen::Vector3d({0, -1, 0}));
    row1.push_back(&n8);
    controlPoints.push_back(row1);

    std::vector<NuTo::DofNode*> row2;
    NuTo::DofNode n9 = NuTo::DofNode(Eigen::Vector3d({0, -1, 1}));
    row2.push_back(&n9);
    NuTo::DofNode n10 = NuTo::DofNode(Eigen::Vector3d({-1, -1, 1}));
    row2.push_back(&n10);
    NuTo::DofNode n11 = NuTo::DofNode(Eigen::Vector3d({-1, 0, 1}));
    row2.push_back(&n11);
    NuTo::DofNode n12 = NuTo::DofNode(Eigen::Vector3d({-1, 1, 1}));
    row2.push_back(&n12);
    NuTo::DofNode n13 = NuTo::DofNode(Eigen::Vector3d({0, 1, 1}));
    row2.push_back(&n13);
    NuTo::DofNode n14 = NuTo::DofNode(Eigen::Vector3d({1, 1, 1}));
    row2.push_back(&n14);
    NuTo::DofNode n15 = NuTo::DofNode(Eigen::Vector3d({1, 0, 1}));
    row2.push_back(&n15);
    NuTo::DofNode n16 = NuTo::DofNode(Eigen::Vector3d({1, -1, 1}));
    row2.push_back(&n16);
    NuTo::DofNode n17 = NuTo::DofNode(Eigen::Vector3d({0, -1, 1}));
    row2.push_back(&n17);
    controlPoints.push_back(row2);

    std::vector<NuTo::DofNode*> row3;
    NuTo::DofNode n18 = NuTo::DofNode(Eigen::Vector3d({0, -1, 2}));
    row3.push_back(&n18);
    NuTo::DofNode n19 = NuTo::DofNode(Eigen::Vector3d({-1, -1, 2}));
    row3.push_back(&n19);
    NuTo::DofNode n20 = NuTo::DofNode(Eigen::Vector3d({-1, 0, 2}));
    row3.push_back(&n20);
    NuTo::DofNode n21 = NuTo::DofNode(Eigen::Vector3d({-1, 1, 2}));
    row3.push_back(&n21);
    NuTo::DofNode n22 = NuTo::DofNode(Eigen::Vector3d({0, 1, 2}));
    row3.push_back(&n22);
    NuTo::DofNode n23 = NuTo::DofNode(Eigen::Vector3d({1, 1, 2}));
    row3.push_back(&n23);
    NuTo::DofNode n24 = NuTo::DofNode(Eigen::Vector3d({1, 0, 2}));
    row3.push_back(&n24);
    NuTo::DofNode n25 = NuTo::DofNode(Eigen::Vector3d({1, -1, 2}));
    row3.push_back(&n25);
    NuTo::DofNode n26 = NuTo::DofNode(Eigen::Vector3d({0, -1, 2}));
    row3.push_back(&n26);
    controlPoints.push_back(row3);

    std::vector<double> knots1 = {0, 0, 0, 1 / 4., 1 / 4., 1 / 2., 1 / 2., 3 / 4., 3 / 4., 1, 1, 1};
    std::vector<double> knots2 = {0, 0, 0, 1, 1, 1};
    std::array<std::vector<double>, 2> knots = {knots1, knots2};

    std::vector<double> weights1 = {1, std::sqrt(2) / 2, 1, std::sqrt(2) / 2, 1, std::sqrt(2) / 2,
                                    1, std::sqrt(2) / 2, 1};
    std::vector<std::vector<double>> weights = {weights1, weights1, weights1};

    std::array<int, 2> degree = {2, 2};

    NuTo::Nurbs<2> surface(knots, controlPoints, weights, degree);

    Eigen::Matrix<double, 2, 1> vec;

    for (int i = 0; i < 100; i++)
    {
        for (int j = 0; j < 100; j++)
        {
            vec << i / 100., j / 100.;
            Eigen::VectorXd vel = surface.Evaluate(vec);
            BOOST_CHECK_CLOSE(vel(0) * vel(0) + vel(1) * vel(1), 1, 1.e-10);
        }
    }

    Eigen::Vector3d derivativeCheck(-1, 0, 0);
    vec << 0, 0;
    Eigen::VectorXd derivative = surface.Evaluate(vec, 1);
    derivative.normalize();
    BoostUnitTest::CheckEigenMatrix(derivative, derivativeCheck);

    vec << 0.25, 0.5;
    derivativeCheck << 0, 1, 0;
    derivative = surface.Evaluate(vec, 1);
    derivative.normalize();
    //    std::cout << derivative << std::endl;
    BoostUnitTest::CheckEigenMatrix(derivative, derivativeCheck);

    vec << 0.5, 0.123;
    derivativeCheck << 1, 0, 0;
    derivative = surface.Evaluate(vec, 1);
    derivative.normalize();
    //    std::cout << derivative << std::endl;
    BoostUnitTest::CheckEigenMatrix(derivative, derivativeCheck);

    vec << 0.75, 0.457;
    derivativeCheck << 0, -1, 0;
    derivative = surface.Evaluate(vec, 1);
    derivative.normalize();
    //    std::cout << derivative << std::endl;
    BoostUnitTest::CheckEigenMatrix(derivative, derivativeCheck);

    vec << 1, 0.891;
    derivativeCheck << -1, 0, 0;
    derivative = surface.Evaluate(vec, 1);
    derivative.normalize();
    BoostUnitTest::CheckEigenMatrix(derivative, derivativeCheck);
}
