#include "BoostUnitTest.h"
#include <boost/test/output_test_stream.hpp>
#include <type_traits>

#include "nuto/mechanics/elements/ElementIga.h"

BOOST_AUTO_TEST_CASE(ElementCopyMove)
{
    BOOST_CHECK(std::is_copy_constructible<NuTo::ElementIga<1>>::value);
    BOOST_CHECK(std::is_move_constructible<NuTo::ElementIga<1>>::value);

    BOOST_CHECK(std::is_copy_constructible<NuTo::ElementIga<2>>::value);
    BOOST_CHECK(std::is_move_constructible<NuTo::ElementIga<2>>::value);
}


BOOST_AUTO_TEST_CASE(ExtractNodeValues1D)
{
    // IGA geometry of a circle (Nurbs curve should exactly fit the circle)
    std::vector<std::vector<NuTo::DofNode*>> controlPoints;

    NuTo::DofNode n1 = NuTo::DofNode(Eigen::Vector2d({0, -1}));
    NuTo::DofNode n2 = NuTo::DofNode(Eigen::Vector2d({-1, -1}));
    NuTo::DofNode n3 = NuTo::DofNode(Eigen::Vector2d({-1, 0}));
    NuTo::DofNode n4 = NuTo::DofNode(Eigen::Vector2d({-1, 1}));
    NuTo::DofNode n5 = NuTo::DofNode(Eigen::Vector2d({0, 1}));
    NuTo::DofNode n6 = NuTo::DofNode(Eigen::Vector2d({1, 1}));
    NuTo::DofNode n7 = NuTo::DofNode(Eigen::Vector2d({1, 0}));
    NuTo::DofNode n8 = NuTo::DofNode(Eigen::Vector2d({1, -1}));
    NuTo::DofNode n9 = NuTo::DofNode(Eigen::Vector2d({0, -1}));
    controlPoints.push_back({&n1, &n2, &n3, &n4, &n5, &n6, &n7, &n8, &n9});

    std::vector<double> knots1D = {0, 0, 0, 1 / 4., 1 / 4., 1 / 2., 1 / 2., 3 / 4., 3 / 4., 1, 1, 1};
    std::array<std::vector<double>, 1> knots = {knots1D};

    std::vector<double> weights1D = {1, std::sqrt(2) / 2, 1, std::sqrt(2) / 2, 1, std::sqrt(2) / 2,
                                     1, std::sqrt(2) / 2, 1};
    std::vector<std::vector<double>> weights = {weights1D};

    std::array<int, 1> degree = {2};

    NuTo::Nurbs<1> curve(knots, controlPoints, weights, degree);

    std::array<int, 1> knotIDsCell = {4};
    NuTo::ElementIga<1> iga(knotIDsCell, curve);
    Eigen::VectorXd nodeValues = iga.ExtractNodeValues();

    // ip coordinates are passed in
    Eigen::VectorXd param(1);
    param << 0.;
    Eigen::VectorXd shapefuns = iga.GetShapeFunctions(param);

    Eigen::Vector2d result;
    result(0) = shapefuns(0) * nodeValues(0) + shapefuns(1) * nodeValues(2) + shapefuns(2) * nodeValues(4);
    result(1) = shapefuns(0) * nodeValues(1) + shapefuns(1) * nodeValues(3) + shapefuns(2) * nodeValues(5);

    // curve parameter coordinates are passed in
    param << 0.375;
    BoostUnitTest::CheckVector(result, curve.Evaluate(param), 2);
}

BOOST_AUTO_TEST_CASE(ExtractNodeValuesInstance)
{
    std::vector<std::vector<NuTo::DofNode*>> controlPoints;
    NuTo::DofNode n1 = NuTo::DofNode(2, 2);
    NuTo::DofNode n2 = NuTo::DofNode(2, 2);

    n1.SetValues(Eigen::Vector2d(11, 12), 1);
    n2.SetValues(Eigen::Vector2d(13, 14), 1);
    controlPoints.push_back({&n1, &n2});

    std::vector<double> knots1D = {0, 0, 1, 1};
    std::vector<double> weights1D = {1, std::sqrt(2) / 2};

    NuTo::Nurbs<1> curve({knots1D}, controlPoints, {weights1D}, {1});

    NuTo::ElementIga<1> iga({1}, curve);
    BoostUnitTest::CheckEigenMatrix(iga.ExtractNodeValues(1), Eigen::Vector4d(11, 12, 13, 14));
}
