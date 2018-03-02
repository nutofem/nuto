#include "BoostUnitTest.h"

#include "mechanics/integrands/NeumannBc.h"

#include "mechanics/nodes/NodeSimple.h"
#include "mechanics/elements/ElementCollection.h"
#include "mechanics/interpolation/InterpolationTrussLinear.h"
#include "mechanics/interpolation/InterpolationTriangleLinear.h"
#include "mechanics/interpolation/InterpolationTriangleQuadratic.h"
#include "mechanics/cell/Matrix.h"
#include "mechanics/elements/ElementShapeFunctions.h"

using namespace NuTo;

constexpr int cellID = 0;
constexpr int ipNum = 0;

BOOST_AUTO_TEST_CASE(NeumannBc1Din2D)
{
    // coordinate element
    NodeSimple nc0(Eigen::Vector2d(0, 0));
    NodeSimple nc1(Eigen::Vector2d(2, 2));
    InterpolationTrussLinear coordinateInterpolation;
    ElementCollectionFem element({{nc0, nc1}, coordinateInterpolation});

    // displacement nodes
    DofType dof("displacements", 2);
    NodeSimple nd0(Eigen::Vector2d(0, 0));
    NodeSimple nd1(Eigen::Vector2d(0, 0));
    InterpolationTrussLinear displacementInterpolation;
    element.AddDofElement(dof, {{nd0, nd1}, displacementInterpolation});

    Eigen::Vector2d p(4, 42);

    Integrands::NeumannBc<2> neumannIntegrand(dof, p);

    CellData cellData(element, cellID);
    Jacobian dummyJac(element.CoordinateElement().ExtractNodeValues(), // jacobian not needed for N.
                      element.CoordinateElement().GetDerivativeShapeFunctions(Eigen::VectorXd::Constant(1, 0.)), 2);


    // Gradient. What should happen here?
    //
    // We are only dealing with N, dof values do not matter.
    // We want a result vector of length 4 containing ( if we speak in forces f )
    // [ f0x, f0y, f1x, f1y ]
    //
    // On integration point (-1) all the force contribution should be at f0 == p
    {
        Eigen::VectorXd ip = Eigen::VectorXd::Constant(1, -1.0);
        CellIpData cellIpData(cellData, dummyJac, ip, ipNum);

        Eigen::Vector4d expected(p[0], p[1], 0, 0);
        auto gradient = neumannIntegrand.ExternalLoad(cellIpData);
        BoostUnitTest::CheckEigenMatrix(gradient[dof], expected);
    }

    // On integration point (+1) all the force contribution should be at f1 == p
    {
        Eigen::VectorXd ip = Eigen::VectorXd::Constant(1, 1.0);
        CellIpData cellIpData(cellData, dummyJac, ip, ipNum);

        Eigen::Vector4d expected(0, 0, p[0], p[1]);
        auto gradient = neumannIntegrand.ExternalLoad(cellIpData);
        BoostUnitTest::CheckEigenMatrix(gradient[dof], expected);
    }

    // On integration point (0) all the force contribution should be spread equally
    {
        Eigen::VectorXd ip = Eigen::VectorXd::Constant(1, 0.0);
        CellIpData cellIpData(cellData, dummyJac, ip, ipNum);

        Eigen::Vector4d expected(p[0] / 2, p[1] / 2, p[0] / 2, p[1] / 2);
        auto gradient = neumannIntegrand.ExternalLoad(cellIpData);
        BoostUnitTest::CheckEigenMatrix(gradient[dof], expected);
    }
}


BOOST_AUTO_TEST_CASE(NeumannBc2Din3D)
{
    std::vector<NodeSimple> nodesDisp;
    for (int i = 0; i < 3; i++)
        nodesDisp.push_back(NodeSimple(Eigen::Vector3d(0, 0, 0)));

    std::vector<NodeSimple> nodes;
    nodes.push_back(NodeSimple(Eigen::Vector3d(1, 0, 0)));
    nodes.push_back(NodeSimple(Eigen::Vector3d(0, 1, 0)));
    nodes.push_back(NodeSimple(Eigen::Vector3d(0, 0, 1)));

    std::vector<NodeSimple*> nodePtrs;
    std::vector<NodeSimple*> nodesDispPtrs;

    nodePtrs.push_back(&nodes[0]);
    nodePtrs.push_back(&nodes[1]);
    nodePtrs.push_back(&nodes[2]);

    nodesDispPtrs.push_back(&nodesDisp[0]);
    nodesDispPtrs.push_back(&nodesDisp[1]);
    nodesDispPtrs.push_back(&nodesDisp[2]);

    // coordinate element
    InterpolationTriangleLinear coordinateInterpolation;
    ElementCollectionFem element({nodePtrs, coordinateInterpolation});

    // displacement nodes
    DofType dof("displacements", 3);
    InterpolationTriangleLinear displacementInterpolation;
    element.AddDofElement(dof, {nodesDispPtrs, displacementInterpolation});

    Eigen::Vector3d p(1., 1., 1.);

    Integrands::NeumannBc<3> neumannIntegrand(dof, p);

    CellData cellData(element, cellID);
    Jacobian dummyJac(element.CoordinateElement().ExtractNodeValues(), // jacobian not needed for N.
                      element.CoordinateElement().GetDerivativeShapeFunctions(Eigen::VectorXd::Constant(2, 1, 0.)), 3);

    std::vector<Eigen::Vector2d> points;
    points.push_back(Eigen::Vector2d(1. / 6., 1. / 6.));
    points.push_back(Eigen::Vector2d(4. / 6., 1. / 6.));
    points.push_back(Eigen::Vector2d(1. / 6., 4. / 6.));

    for (int ip = 0; ip < 3; ip++)
    {
        CellIpData cellIpData(cellData, dummyJac, points[ip], ipNum);

        auto gradient = neumannIntegrand.ExternalLoad(cellIpData);

        Eigen::VectorXd referenceResultOnIP = (1. / 6.) * Eigen::VectorXd::Ones(9);
        referenceResultOnIP.block(3 * ip, 0, 3, 1) = (2. / 3.) * Eigen::VectorXd::Ones(3);

        // results on each ip
        BoostUnitTest::CheckEigenMatrix(gradient[dof], referenceResultOnIP);
    }
}
