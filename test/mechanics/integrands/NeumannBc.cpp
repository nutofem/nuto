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

BOOST_AUTO_TEST_CASE(NeumannBc1Din2D)
{
    // coordinate element
    NodeSimple nc0(Eigen::Vector2d(0, 0));
    NodeSimple nc1(Eigen::Vector2d(2, 2));
    InterpolationTrussLinear coordinateInterpolation(2);
    ElementCollectionFem element({{nc0, nc1}, coordinateInterpolation});

    // displacement nodes
    DofType dof("displacements", 2);
    NodeSimple nd0(Eigen::Vector2d(0, 0));
    NodeSimple nd1(Eigen::Vector2d(0, 0));
    InterpolationTrussLinear displacementInterpolation(2);
    element.AddDofElement(dof, {{nd0, nd1}, displacementInterpolation});

    Eigen::Vector2d p(4, 42);

    Integrands::TimeDependent::NeumannBc<2> neumannIntegrand(dof, p);

    CellData cellData(element);
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
        CellIpData cellIpData(element, dummyJac, ip);

        Eigen::Vector4d expected(p[0], p[1], 0, 0);
        auto gradient = neumannIntegrand.Gradient(cellData, cellIpData);
        BoostUnitTest::CheckEigenMatrix(gradient[dof], expected);
    }

    // On integration point (+1) all the force contribution should be at f1 == p
    {
        Eigen::VectorXd ip = Eigen::VectorXd::Constant(1, 1.0);
        CellIpData cellIpData(element, dummyJac, ip);

        Eigen::Vector4d expected(0, 0, p[0], p[1]);
        auto gradient = neumannIntegrand.Gradient(cellData, cellIpData);
        BoostUnitTest::CheckEigenMatrix(gradient[dof], expected);
    }

    // On integration point (0) all the force contribution should be spread equally
    {
        Eigen::VectorXd ip = Eigen::VectorXd::Constant(1, 0.0);
        CellIpData cellIpData(element, dummyJac, ip);

        Eigen::Vector4d expected(p[0] / 2, p[1] / 2, p[0] / 2, p[1] / 2);
        auto gradient = neumannIntegrand.Gradient(cellData, cellIpData);
        BoostUnitTest::CheckEigenMatrix(gradient[dof], expected);
    }

    // Hessian. What should happen here?
    // Zeros in the right dimension

    CellIpData cellIpData(element, dummyJac, Eigen::VectorXd::Constant(1, 0));
    auto hessian0 = neumannIntegrand.Hessian0(cellData, cellIpData);
    BoostUnitTest::CheckEigenMatrix(hessian0(dof, dof), Eigen::MatrixXd::Zero(4, 4));
}

BOOST_AUTO_TEST_CASE(NeumannBc2Din3D_TriangleLinear)
{
    std::vector<NodeSimple> nodesDisp;
    for (int i = 0; i < 4; i++)
        nodesDisp.push_back(NodeSimple(Eigen::Vector3d(0, 0, 0)));

    std::vector<NodeSimple> nodes;
    nodes.push_back(NodeSimple(Eigen::Vector3d(1, 0, 0)));
    nodes.push_back(NodeSimple(Eigen::Vector3d(0, 1, 0)));
    nodes.push_back(NodeSimple(Eigen::Vector3d(0, 0, 1)));
    nodes.push_back(NodeSimple(Eigen::Vector3d(1, 1, 1)));

    std::vector<NodeSimple*> nodePtrs;
    std::vector<NodeSimple*> nodesDispPtrs;

    std::vector<Eigen::Vector2d> points;
    points.push_back(Eigen::Vector2d(4. / 6., 1. / 6.));
    points.push_back(Eigen::Vector2d(1. / 6., 1. / 6.));
    points.push_back(Eigen::Vector2d(1. / 6., 4. / 6.));

    for (Eigen::Vector2d& surfaceIP : points)
    {
        for (int surface = 0; surface < 4; surface++)
        {
            Eigen::Matrix<double, 9, 12> transform3d2d = Eigen::Matrix<double, 9, 12>::Zero();
            Eigen::Vector3d volumetricParameters;
            nodePtrs.clear();
            nodesDispPtrs.clear();
            switch (surface)
            {
            case 0:
            {
                volumetricParameters << surfaceIP(0), surfaceIP(1), 0; // x,y
                nodePtrs.push_back(&nodes[0]);
                nodePtrs.push_back(&nodes[1]);
                nodePtrs.push_back(&nodes[2]);

                nodesDispPtrs.push_back(&nodesDisp[0]);
                nodesDispPtrs.push_back(&nodesDisp[1]);
                nodesDispPtrs.push_back(&nodesDisp[2]);

                transform3d2d.block(0, 0, 3, 3) = Eigen::Matrix<double, 3, 3>::Identity();
                transform3d2d.block(3, 3, 3, 3) = Eigen::Matrix<double, 3, 3>::Identity();
                transform3d2d.block(6, 6, 3, 3) = Eigen::Matrix<double, 3, 3>::Identity();

                break;
            }
            case 1:
            {
                volumetricParameters << 0, surfaceIP(0), surfaceIP(1); // y, z
                nodePtrs.push_back(&nodes[0]);
                nodePtrs.push_back(&nodes[2]);
                nodePtrs.push_back(&nodes[3]);

                nodesDispPtrs.push_back(&nodesDisp[0]);
                nodesDispPtrs.push_back(&nodesDisp[2]);
                nodesDispPtrs.push_back(&nodesDisp[3]);

                transform3d2d.block(0, 0, 3, 3) = Eigen::Matrix<double, 3, 3>::Identity();
                transform3d2d.block(3, 6, 3, 3) = Eigen::Matrix<double, 3, 3>::Identity();
                transform3d2d.block(6, 9, 3, 3) = Eigen::Matrix<double, 3, 3>::Identity();
                break;
            }
            case 2:
            {
                volumetricParameters << surfaceIP(0), 0, surfaceIP(1); // x, z
                nodePtrs.push_back(&nodes[0]);
                nodePtrs.push_back(&nodes[1]);
                nodePtrs.push_back(&nodes[3]);

                nodesDispPtrs.push_back(&nodesDisp[0]);
                nodesDispPtrs.push_back(&nodesDisp[1]);
                nodesDispPtrs.push_back(&nodesDisp[3]);

                transform3d2d.block(0, 0, 3, 3) = Eigen::Matrix<double, 3, 3>::Identity();
                transform3d2d.block(3, 3, 3, 3) = Eigen::Matrix<double, 3, 3>::Identity();
                transform3d2d.block(6, 9, 3, 3) = Eigen::Matrix<double, 3, 3>::Identity();
                break;
            }
            case 3:
            {
                volumetricParameters << 1 - surfaceIP(0) - surfaceIP(1), surfaceIP(0), surfaceIP(1); // x, y, z
                nodePtrs.push_back(&nodes[1]);
                nodePtrs.push_back(&nodes[2]);
                nodePtrs.push_back(&nodes[3]);

                nodesDispPtrs.push_back(&nodesDisp[1]);
                nodesDispPtrs.push_back(&nodesDisp[2]);
                nodesDispPtrs.push_back(&nodesDisp[3]);

                transform3d2d.block(0, 3, 3, 3) = Eigen::Matrix<double, 3, 3>::Identity();
                transform3d2d.block(3, 6, 3, 3) = Eigen::Matrix<double, 3, 3>::Identity();
                transform3d2d.block(6, 9, 3, 3) = Eigen::Matrix<double, 3, 3>::Identity();
                break;
            }
            }

            // coordinate element
            InterpolationTriangleLinear coordinateInterpolation(3);
            ElementCollectionFem element({nodePtrs, coordinateInterpolation});

            // displacement nodes
            DofType dof("displacements", 3);
            InterpolationTriangleLinear displacementInterpolation(3);
            element.AddDofElement(dof, {nodesDispPtrs, displacementInterpolation});

            Eigen::Vector3d p(1, 1, 2);

            Integrands::TimeDependent::NeumannBc<3> neumannIntegrand(dof, p);

            CellData cellData(element);
            Jacobian dummyJac(
                    element.CoordinateElement().ExtractNodeValues(), // jacobian not needed for N.
                    element.CoordinateElement().GetDerivativeShapeFunctions(Eigen::VectorXd::Constant(2, 1, 0.)), 3);


            CellIpData cellIpData(element, dummyJac, surfaceIP);

            Eigen::Matrix<double, 4, 1> TetrahedronDN =
                    NuTo::ShapeFunctions3D::ShapeFunctionsTetrahedronOrder1(volumetricParameters);

            Eigen::MatrixXd NTetrahedron = NuTo::Matrix::N(TetrahedronDN, 4, 3);
            auto gradient = neumannIntegrand.Gradient(cellData, cellIpData);
            BoostUnitTest::CheckEigenMatrix(gradient[dof], transform3d2d * NTetrahedron.transpose() * p);
        }
    }
}


BOOST_AUTO_TEST_CASE(NeumannBc2Din3D_TriangleQuadratic)
{
    std::vector<NodeSimple> nodesDisp;
    for (int i = 0; i < 10; i++)
        nodesDisp.push_back(NodeSimple(Eigen::Vector3d(0, 0, 0)));

    std::vector<NodeSimple> nodes;
    nodes.push_back(NodeSimple(Eigen::Vector3d(1, 0, 0)));
    nodes.push_back(NodeSimple(Eigen::Vector3d(0, 1, 0)));
    nodes.push_back(NodeSimple(Eigen::Vector3d(0, 0, 1)));
    nodes.push_back(NodeSimple(Eigen::Vector3d(1, 1, 1)));
    nodes.push_back(NodeSimple(Eigen::Vector3d(0.5, 0.5, 0)));
    nodes.push_back(NodeSimple(Eigen::Vector3d(0, 0.5, 0.5)));
    nodes.push_back(NodeSimple(Eigen::Vector3d(0.5, 0, 0.5)));
    nodes.push_back(NodeSimple(Eigen::Vector3d(1, 0.5, 0.5)));
    nodes.push_back(NodeSimple(Eigen::Vector3d(0.5, 0.5, 1)));
    nodes.push_back(NodeSimple(Eigen::Vector3d(0.5, 1, 0.5)));


    std::vector<NodeSimple*> nodePtrs;
    std::vector<NodeSimple*> nodesDispPtrs;

    for (NodeSimple it : nodes)
        nodePtrs.push_back(&it);

    for (NodeSimple it : nodesDisp)
        nodesDispPtrs.push_back(&it);

    std::vector<Eigen::Vector2d> points;
    points.push_back(Eigen::Vector2d(1. / 6., 1. / 6.));
    points.push_back(Eigen::Vector2d(4. / 6., 1. / 6.));
    points.push_back(Eigen::Vector2d(1. / 6., 4. / 6.));


    for (Eigen::Vector2d& surfaceIP : points)
    {

        for (int surface = 0; surface < 4; surface++)
        {
            Eigen::Vector3d volumetricParameters;
            Eigen::Matrix<double, 18, 30> transform3d2d = Eigen::Matrix<double, 18, 30>::Zero();
            nodePtrs.clear();
            nodesDispPtrs.clear();

            switch (surface)
            {
            case 0:
            {
                volumetricParameters << surfaceIP(0), surfaceIP(1), 0; // x,y
                nodePtrs.push_back(&nodes[0]);
                nodePtrs.push_back(&nodes[1]);
                nodePtrs.push_back(&nodes[2]);

                nodePtrs.push_back(&nodes[4]);
                nodePtrs.push_back(&nodes[5]);
                nodePtrs.push_back(&nodes[6]);

                nodesDispPtrs.push_back(&nodesDisp[0]);
                nodesDispPtrs.push_back(&nodesDisp[1]);
                nodesDispPtrs.push_back(&nodesDisp[2]);

                nodesDispPtrs.push_back(&nodesDisp[4]);
                nodesDispPtrs.push_back(&nodesDisp[5]);
                nodesDispPtrs.push_back(&nodesDisp[6]);

                transform3d2d.block(0 * 3, 0 * 3, 3, 3) = Eigen::Matrix<double, 3, 3>::Identity();
                transform3d2d.block(1 * 3, 1 * 3, 3, 3) = Eigen::Matrix<double, 3, 3>::Identity();
                transform3d2d.block(2 * 3, 2 * 3, 3, 3) = Eigen::Matrix<double, 3, 3>::Identity();

                transform3d2d.block(3 * 3, 4 * 3, 3, 3) = Eigen::Matrix<double, 3, 3>::Identity();
                transform3d2d.block(4 * 3, 5 * 3, 3, 3) = Eigen::Matrix<double, 3, 3>::Identity();
                transform3d2d.block(5 * 3, 6 * 3, 3, 3) = Eigen::Matrix<double, 3, 3>::Identity();


                break;
            }
            case 1:
            {
                volumetricParameters << 0, surfaceIP(0), surfaceIP(1); // y, z
                nodePtrs.push_back(&nodes[0]);
                nodePtrs.push_back(&nodes[2]);
                nodePtrs.push_back(&nodes[3]);

                nodePtrs.push_back(&nodes[6]);
                nodePtrs.push_back(&nodes[8]);
                nodePtrs.push_back(&nodes[7]);

                nodesDispPtrs.push_back(&nodesDisp[0]);
                nodesDispPtrs.push_back(&nodesDisp[2]);
                nodesDispPtrs.push_back(&nodesDisp[3]);

                nodesDispPtrs.push_back(&nodesDisp[6]);
                nodesDispPtrs.push_back(&nodesDisp[8]);
                nodesDispPtrs.push_back(&nodesDisp[7]);


                transform3d2d.block(0 * 3, 0 * 3, 3, 3) = Eigen::Matrix<double, 3, 3>::Identity();
                transform3d2d.block(1 * 3, 2 * 3, 3, 3) = Eigen::Matrix<double, 3, 3>::Identity();
                transform3d2d.block(2 * 3, 3 * 3, 3, 3) = Eigen::Matrix<double, 3, 3>::Identity();

                transform3d2d.block(3 * 3, 6 * 3, 3, 3) = Eigen::Matrix<double, 3, 3>::Identity();
                transform3d2d.block(4 * 3, 8 * 3, 3, 3) = Eigen::Matrix<double, 3, 3>::Identity();
                transform3d2d.block(5 * 3, 7 * 3, 3, 3) = Eigen::Matrix<double, 3, 3>::Identity();
                break;
            }
            case 2:
            {
                volumetricParameters << surfaceIP(0), 0, surfaceIP(1); // x, z
                nodePtrs.push_back(&nodes[0]);
                nodePtrs.push_back(&nodes[1]);
                nodePtrs.push_back(&nodes[3]);

                nodePtrs.push_back(&nodes[4]);
                nodePtrs.push_back(&nodes[9]);
                nodePtrs.push_back(&nodes[7]);

                nodesDispPtrs.push_back(&nodesDisp[0]);
                nodesDispPtrs.push_back(&nodesDisp[1]);
                nodesDispPtrs.push_back(&nodesDisp[3]);

                nodesDispPtrs.push_back(&nodesDisp[4]);
                nodesDispPtrs.push_back(&nodesDisp[9]);
                nodesDispPtrs.push_back(&nodesDisp[7]);

                transform3d2d.block(0 * 3, 0 * 3, 3, 3) = Eigen::Matrix<double, 3, 3>::Identity();
                transform3d2d.block(1 * 3, 1 * 3, 3, 3) = Eigen::Matrix<double, 3, 3>::Identity();
                transform3d2d.block(2 * 3, 3 * 3, 3, 3) = Eigen::Matrix<double, 3, 3>::Identity();

                transform3d2d.block(3 * 3, 4 * 3, 3, 3) = Eigen::Matrix<double, 3, 3>::Identity();
                transform3d2d.block(4 * 3, 9 * 3, 3, 3) = Eigen::Matrix<double, 3, 3>::Identity();
                transform3d2d.block(5 * 3, 7 * 3, 3, 3) = Eigen::Matrix<double, 3, 3>::Identity();
                break;
            }
            case 3:
            {
                volumetricParameters << 1 - surfaceIP(0) - surfaceIP(1), surfaceIP(0), surfaceIP(1); // x, y, z
                nodePtrs.push_back(&nodes[1]);
                nodePtrs.push_back(&nodes[2]);
                nodePtrs.push_back(&nodes[3]);

                nodePtrs.push_back(&nodes[5]);
                nodePtrs.push_back(&nodes[8]);
                nodePtrs.push_back(&nodes[9]);

                nodesDispPtrs.push_back(&nodesDisp[1]);
                nodesDispPtrs.push_back(&nodesDisp[2]);
                nodesDispPtrs.push_back(&nodesDisp[3]);

                nodesDispPtrs.push_back(&nodesDisp[5]);
                nodesDispPtrs.push_back(&nodesDisp[8]);
                nodesDispPtrs.push_back(&nodesDisp[9]);

                transform3d2d.block(0 * 3, 1 * 3, 3, 3) = Eigen::Matrix<double, 3, 3>::Identity();
                transform3d2d.block(1 * 3, 2 * 3, 3, 3) = Eigen::Matrix<double, 3, 3>::Identity();
                transform3d2d.block(2 * 3, 3 * 3, 3, 3) = Eigen::Matrix<double, 3, 3>::Identity();

                transform3d2d.block(3 * 3, 5 * 3, 3, 3) = Eigen::Matrix<double, 3, 3>::Identity();
                transform3d2d.block(4 * 3, 8 * 3, 3, 3) = Eigen::Matrix<double, 3, 3>::Identity();
                transform3d2d.block(5 * 3, 9 * 3, 3, 3) = Eigen::Matrix<double, 3, 3>::Identity();
                break;
            }
            }

            // coordinate element
            InterpolationTriangleQuadratic coordinateInterpolation(3);
            ElementCollectionFem element({nodePtrs, coordinateInterpolation});

            // displacement nodes
            DofType dof("displacements", 3);
            InterpolationTriangleQuadratic displacementInterpolation(3);
            element.AddDofElement(dof, {nodesDispPtrs, displacementInterpolation});

            Eigen::Vector3d p(1, 1, 1);

            Integrands::TimeDependent::NeumannBc<3> neumannIntegrand(dof, p);

            CellData cellData(element);
            Jacobian dummyJac(
                    element.CoordinateElement().ExtractNodeValues(), // jacobian not needed for N.
                    element.CoordinateElement().GetDerivativeShapeFunctions(Eigen::VectorXd::Constant(2, 1, 0.)), 3);

            CellIpData cellIpData(element, dummyJac, surfaceIP);

            Eigen::VectorXd TetrahedronDN =
                    NuTo::ShapeFunctions3D::ShapeFunctionsTetrahedronOrder2(volumetricParameters);

            Eigen::MatrixXd NTetrahedron = NuTo::Matrix::N(TetrahedronDN, 10, 3);
            auto gradient = neumannIntegrand.Gradient(cellData, cellIpData);
            BoostUnitTest::CheckEigenMatrix(gradient[dof], transform3d2d * NTetrahedron.transpose() * p);
        }
    }
}
