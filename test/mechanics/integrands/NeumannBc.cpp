#include "BoostUnitTest.h"

#include "mechanics/integrands/NeumannBc.h"

#include "mechanics/nodes/NodeSimple.h"
#include "mechanics/elements/ElementCollection.h"
#include "mechanics/interpolation/InterpolationTrussLinear.h"

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

    Integrands::NeumannBc<2> neumannIntegrand(dof, p);

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
