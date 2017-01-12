/*
 * Interpolation3DPrism.cpp
 *
 *  Created on: 10 January 2017
 *      Author: Thomas Titscher
 */

#include "mechanics/MechanicsException.h"
#include "mechanics/elements/ElementShapeFunctions.h"
#include "mechanics/integrationtypes/IntegrationTypeEnum.h"
#include "mechanics/interpolationtypes/InterpolationTypeEnum.h"
#include "mechanics/interpolationtypes/Interpolation3DPrism.h"

NuTo::Interpolation3DPrism::Interpolation3DPrism(NuTo::Node::eDof rDofType, NuTo::Interpolation::eTypeOrder rTypeOrder, int rDimension) :
        Interpolation3D::Interpolation3D(rDofType, rTypeOrder, rDimension)
{
    Initialize();
}

NuTo::eIntegrationType NuTo::Interpolation3DPrism::GetStandardIntegrationType() const
{
    switch (mTypeOrder)
    {
    case NuTo::Interpolation::eTypeOrder::EQUIDISTANT1:
        return NuTo::eIntegrationType::IntegrationType3D6NGauss2x3Ip;
    case NuTo::Interpolation::eTypeOrder::EQUIDISTANT2:
        return NuTo::eIntegrationType::IntegrationType3D6NGauss2x3Ip;
    default:
        throw MechanicsException(__PRETTY_FUNCTION__, "Interpolation for exact integration of " + Interpolation::TypeOrderToString(mTypeOrder) + " not implemented");
    }
}

Eigen::VectorXd NuTo::Interpolation3DPrism::CalculateShapeFunctions(const Eigen::VectorXd& rCoordinates) const
{
    switch (mTypeOrder)
    {
    case NuTo::Interpolation::eTypeOrder::EQUIDISTANT1:
        return ShapeFunctions3D::ShapeFunctionsPrismOrder1(rCoordinates);
    case NuTo::Interpolation::eTypeOrder::EQUIDISTANT2:
        return ShapeFunctions3D::ShapeFunctionsPrismOrder2(rCoordinates);
    default:
        throw MechanicsException(__PRETTY_FUNCTION__, "Interpolation order for " + Interpolation::TypeOrderToString(mTypeOrder) + " not implemented");
    }
}

Eigen::MatrixXd NuTo::Interpolation3DPrism::CalculateDerivativeShapeFunctionsNatural(const Eigen::VectorXd& rCoordinates) const
{
    switch (mTypeOrder)
    {
    case NuTo::Interpolation::eTypeOrder::EQUIDISTANT1:
        return ShapeFunctions3D::DerivativeShapeFunctionsPrismOrder1(rCoordinates);
    case NuTo::Interpolation::eTypeOrder::EQUIDISTANT2:
        return ShapeFunctions3D::DerivativeShapeFunctionsPrismOrder2(rCoordinates);
    default:
        throw MechanicsException(__PRETTY_FUNCTION__, "Interpolation order for " + Interpolation::TypeOrderToString(mTypeOrder) + " not implemented");
    }
}

Eigen::VectorXd NuTo::Interpolation3DPrism::CalculateNaturalNodeCoordinates(int rNodeIndexDof) const
{
    switch (mTypeOrder)
    {
    case NuTo::Interpolation::eTypeOrder::EQUIDISTANT1:
        return ShapeFunctions3D::NodeCoordinatesPrismOrder1(rNodeIndexDof);
    case NuTo::Interpolation::eTypeOrder::EQUIDISTANT2:
        return ShapeFunctions3D::NodeCoordinatesPrismOrder2(rNodeIndexDof);
    default:
        throw MechanicsException(__PRETTY_FUNCTION__, "Node arrangement for " + Interpolation::TypeOrderToString(mTypeOrder) + " not implemented");
    }
}

int NuTo::Interpolation3DPrism::CalculateNumNodes() const
{
    switch (mTypeOrder)
    {
    case NuTo::Interpolation::eTypeOrder::EQUIDISTANT1:
        return 6;
    case NuTo::Interpolation::eTypeOrder::EQUIDISTANT2:
        return 15;
    default:
        throw MechanicsException(__PRETTY_FUNCTION__, "Interpolation type and order " + Interpolation::TypeOrderToString(mTypeOrder) + " not implemented");
    }
}

Eigen::VectorXd NuTo::Interpolation3DPrism::CalculateNaturalSurfaceCoordinates(const Eigen::VectorXd& rNaturalSurfaceCoordinates, int rSurface) const
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
    default:
        throw MechanicsException(__PRETTY_FUNCTION__, "BRICK3D has exactly six surfaces, 0 to 5. You tried to access " + std::to_string(rSurface) + ".");
    }
}

Eigen::MatrixXd NuTo::Interpolation3DPrism::CalculateDerivativeNaturalSurfaceCoordinates(const Eigen::VectorXd& rNaturalSurfaceCoordinates, int rSurface) const
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

    default:
        throw MechanicsException(__PRETTY_FUNCTION__, "PRISM3D has exactly five surfaces, 0 to 4. You tried to access " + std::to_string(rSurface) + ".");
    }
    return dXidAlpha;
}

std::vector<Eigen::VectorXd> NuTo::Interpolation3DPrism::GetSurfaceEdgesCoordinates(int rSurface) const
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
