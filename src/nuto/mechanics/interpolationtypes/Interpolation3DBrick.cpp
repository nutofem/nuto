/*
 * Interpolation3DBrick.h
 *
 *  Created on: 19 May 2015
 *      Author: ttitsche
 */

#include "nuto/mechanics/interpolationtypes/Interpolation3DBrick.h"

NuTo::Interpolation3DBrick::Interpolation3DBrick(NuTo::Node::eAttributes rDofType, NuTo::Interpolation::eTypeOrder rTypeOrder) :
        Interpolation3D::Interpolation3D(rDofType, rTypeOrder)
{
    Initialize();
}

NuTo::IntegrationType::eIntegrationType NuTo::Interpolation3DBrick::GetStandardIntegrationType() const
{
    switch (mTypeOrder)
    {
    case NuTo::Interpolation::eTypeOrder::EQUIDISTANT1:
        return NuTo::IntegrationType::IntegrationType3D8NGauss2x2x2Ip;
    case NuTo::Interpolation::eTypeOrder::EQUIDISTANT2:
        return NuTo::IntegrationType::IntegrationType3D8NGauss2x2x2Ip;
    case NuTo::Interpolation::eTypeOrder::LOBATTO2:
        return NuTo::IntegrationType::IntegrationType3D8NLobatto3x3x3Ip;
    case NuTo::Interpolation::eTypeOrder::LOBATTO3:
        return NuTo::IntegrationType::IntegrationType3D8NLobatto4x4x4Ip;
    case NuTo::Interpolation::eTypeOrder::LOBATTO4:
        return NuTo::IntegrationType::IntegrationType3D8NLobatto5x5x5Ip;
    default:
        throw MechanicsException("[NuTo::Interpolation3DBrick::GetStandardIntegrationType] Interpolation for exact integration of " + Interpolation::TypeOrderToString(mTypeOrder) + " not implemented");
    }
}

const Eigen::VectorXd NuTo::Interpolation3DBrick::CalculateShapeFunctions(const Eigen::VectorXd& rCoordinates) const
{
    switch (mTypeOrder)
    {
    case NuTo::Interpolation::eTypeOrder::EQUIDISTANT1:
        return ShapeFunctions3D::ShapeFunctionsBrickOrder1(rCoordinates);
    case NuTo::Interpolation::eTypeOrder::EQUIDISTANT2:
        return ShapeFunctions3D::ShapeFunctionsBrickOrder2(rCoordinates);
    case NuTo::Interpolation::eTypeOrder::LOBATTO2:
        return ShapeFunctions3D::ShapeFunctionsBrickSpectralOrder2(rCoordinates);
    case NuTo::Interpolation::eTypeOrder::LOBATTO3:
        return ShapeFunctions3D::ShapeFunctionsBrickSpectralOrder3(rCoordinates);
    case NuTo::Interpolation::eTypeOrder::LOBATTO4:
        return ShapeFunctions3D::ShapeFunctionsBrickSpectralOrder4(rCoordinates);
    default:
        throw MechanicsException("[NuTo::Interpolation3DBrick::CalculateShapeFunctions] Interpolation order for " + Interpolation::TypeOrderToString(mTypeOrder) + " not implemented");
    }
}

const Eigen::MatrixXd NuTo::Interpolation3DBrick::CalculateDerivativeShapeFunctionsNatural(const Eigen::VectorXd& rCoordinates) const
{
    switch (mTypeOrder)
    {
    case NuTo::Interpolation::eTypeOrder::EQUIDISTANT1:
        return ShapeFunctions3D::DerivativeShapeFunctionsBrickOrder1(rCoordinates);
    case NuTo::Interpolation::eTypeOrder::EQUIDISTANT2:
        return ShapeFunctions3D::DerivativeShapeFunctionsBrickOrder2(rCoordinates);
    case NuTo::Interpolation::eTypeOrder::LOBATTO2:
        return ShapeFunctions3D::DerivativeShapeFunctionsBrickSpectralOrder2(rCoordinates);
    case NuTo::Interpolation::eTypeOrder::LOBATTO3:
        return ShapeFunctions3D::DerivativeShapeFunctionsBrickSpectralOrder3(rCoordinates);
    case NuTo::Interpolation::eTypeOrder::LOBATTO4:
        return ShapeFunctions3D::DerivativeShapeFunctionsBrickSpectralOrder4(rCoordinates);
    default:
        throw MechanicsException("[NuTo::Interpolation3DBrick::CalculateDerivativeShapeFunctionsNatural] Interpolation order for " + Interpolation::TypeOrderToString(mTypeOrder) + " not implemented");
    }
}

const Eigen::VectorXd NuTo::Interpolation3DBrick::CalculateNaturalNodeCoordinates(int rNodeIndexDof) const
{
    switch (mTypeOrder)
    {
    case NuTo::Interpolation::eTypeOrder::EQUIDISTANT1:
        return ShapeFunctions3D::NodeCoordinatesBrickOrder1(rNodeIndexDof);
    case NuTo::Interpolation::eTypeOrder::EQUIDISTANT2:
        return ShapeFunctions3D::NodeCoordinatesBrickOrder2(rNodeIndexDof);
    case NuTo::Interpolation::eTypeOrder::LOBATTO2:
        return ShapeFunctions3D::NodeCoordinatesBrickSpectralOrder2(rNodeIndexDof);
    case NuTo::Interpolation::eTypeOrder::LOBATTO3:
        return ShapeFunctions3D::NodeCoordinatesBrickSpectralOrder3(rNodeIndexDof);
    case NuTo::Interpolation::eTypeOrder::LOBATTO4:
        return ShapeFunctions3D::NodeCoordinatesBrickSpectralOrder4(rNodeIndexDof);
    default:
        throw MechanicsException("[NuTo::Interpolation3DBrick::CalculateNaturalNodeCoordinates] Node arrangement for " + Interpolation::TypeOrderToString(mTypeOrder) + " not implemented");
    }
}

int NuTo::Interpolation3DBrick::CalculateNumNodes() const
{
    switch (mTypeOrder)
    {
    case NuTo::Interpolation::eTypeOrder::EQUIDISTANT1:
        return 8;
    case NuTo::Interpolation::eTypeOrder::EQUIDISTANT2:
        return 20;
    case NuTo::Interpolation::eTypeOrder::LOBATTO2:
        return 27;
    case NuTo::Interpolation::eTypeOrder::LOBATTO3:
        return 64;
    case NuTo::Interpolation::eTypeOrder::LOBATTO4:
        return 125;
    default:
        throw MechanicsException("[NuTo::Interpolation3DBrick::CalculateNumNodes] Interpolation type and order " + Interpolation::TypeOrderToString(mTypeOrder) + " not implemented");
    }
}

const Eigen::VectorXd NuTo::Interpolation3DBrick::CalculateNaturalSurfaceCoordinates(const Eigen::VectorXd& rNaturalSurfaceCoordinates, int rSurface) const
{
    assert(rNaturalSurfaceCoordinates.rows() == 2);
    double alpha = rNaturalSurfaceCoordinates(0);
    double beta = rNaturalSurfaceCoordinates(1);

    switch (rSurface)
    {
    case 0:
        return Eigen::Vector3d( beta,alpha,  -1.);
    case 1:
        return Eigen::Vector3d(alpha,  -1., beta);
    case 2:
        return Eigen::Vector3d(  -1., beta,alpha);
    case 3:
        return Eigen::Vector3d(alpha, beta,   1.);
    case 4:
        return Eigen::Vector3d( beta,   1.,alpha);
    case 5:
        return Eigen::Vector3d(   1.,alpha, beta);
    default:
        throw MechanicsException("[NuTo::Interpolation3DBrick::CalculateNaturalSurfaceCoordinates] BRICK3D has exactly six surfaces, 0 to 5. You tried to access " + std::to_string(rSurface) + ".");
    }
}

const Eigen::MatrixXd NuTo::Interpolation3DBrick::CalculateDerivativeNaturalSurfaceCoordinates(const Eigen::VectorXd& rNaturalSurfaceCoordinates, int rSurface) const
{
    assert(rNaturalSurfaceCoordinates.rows() == 2);
    Eigen::MatrixXd dXidAlpha = Eigen::Matrix<double, 3, 2>::Zero();
    switch (rSurface)
    {
    case 0:
        dXidAlpha(0, 1) = 1.;
        dXidAlpha(1, 0) = 1.;
        break;
    case 1:
        dXidAlpha(0, 0) = 1.;
        dXidAlpha(2, 1) = 1.;
        break;
    case 2:
        dXidAlpha(1, 1) = 1.;
        dXidAlpha(2, 0) = 1.;
        break;
    case 3:
        dXidAlpha(0, 0) = 1.;
        dXidAlpha(1, 1) = 1.;
        break;
    case 4:
        dXidAlpha(0, 1) = 1.;
        dXidAlpha(2, 0) = 1.;
        break;
    case 5:
        dXidAlpha(1, 0) = 1.;
        dXidAlpha(2, 1) = 1.;
        break;

    default:
        throw MechanicsException("[NuTo::Interpolation3DBrick::CalculateDerivativeNaturalSurfaceCoordinates] BRICK3D has exactly six surfaces, 0 to 5. You tried to access " + std::to_string(rSurface) + ".");
    }
    return dXidAlpha;
}

const std::vector<Eigen::VectorXd> NuTo::Interpolation3DBrick::GetSurfaceEdgesCoordinates(int rSurface) const
{
    int numNodes = 4;
    // returns exactly three nodes, (-1,-1).T; (1,-1).T; (1,1).T and (-1,1).T
    Eigen::Matrix<double, 2, 4> alpha = Eigen::Matrix<double, 2, 4>::Zero(); // row1 = alpha, row2 = beta

    alpha << -1.,  1.,  1., -1.,
             -1., -1.,  1.,  1.;
    std::vector<Eigen::VectorXd> surfaceEdgeCoordinates(numNodes);
    for (int i = 0; i < numNodes; ++i)
    {
        Eigen::VectorXd naturalSurfaceCoordinate = alpha.col(i);
        surfaceEdgeCoordinates[i] = CalculateNaturalSurfaceCoordinates(naturalSurfaceCoordinate, rSurface);
    }
    return surfaceEdgeCoordinates;
}
