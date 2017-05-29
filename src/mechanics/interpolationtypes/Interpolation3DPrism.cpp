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
        return NuTo::eIntegrationType::IntegrationType3D6NGauss1Ip;
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
        return 18;
    default:
        throw MechanicsException(__PRETTY_FUNCTION__, "Interpolation type and order " + Interpolation::TypeOrderToString(mTypeOrder) + " not implemented");
    }
}

Eigen::VectorXd NuTo::Interpolation3DPrism::CalculateNaturalSurfaceCoordinates(const Eigen::VectorXd& rNaturalSurfaceCoordinates, int rSurface) const
{
    assert(rNaturalSurfaceCoordinates.rows() == 2);
    double alpha = rNaturalSurfaceCoordinates[0];
    double beta = rNaturalSurfaceCoordinates[1];

    switch (rSurface)
    {
    case 0:
        return Eigen::Vector3d( beta,alpha,  -1.);
    case 1:
        return Eigen::Vector3d(alpha, beta, 1.);
    case 2:
        return Eigen::Vector3d(0.5 + 0.5 * alpha, 0, beta);
    case 3:
        return Eigen::Vector3d(0, 0.5 + 0.5*beta,   alpha);
    case 4:
        return Eigen::Vector3d(0.5 + 0.5*beta, 0.5 - 0.5*beta, alpha);
    default:
        throw MechanicsException(__PRETTY_FUNCTION__, "PRISM3D has exactly five surfaces, 0 to 4. You tried to access " + std::to_string(rSurface) + ".");
    }
}

Eigen::MatrixXd NuTo::Interpolation3DPrism::CalculateDerivativeNaturalSurfaceCoordinates(const Eigen::VectorXd& rNaturalSurfaceCoordinates, int rSurface) const
{
    assert(rNaturalSurfaceCoordinates.rows() == 2);
    Eigen::MatrixXd dXidAlpha = Eigen::Matrix<double, 3, 2>::Zero();
    switch (rSurface)
    {
    case 0:
        dXidAlpha(1, 0) = 1.;
        dXidAlpha(0, 1) = 1.;
        break;
    case 1:
        dXidAlpha(0, 0) = 1.;
        dXidAlpha(1, 1) = 1.;
        break;
    case 2:
        dXidAlpha(0, 0) = 0.5;
        dXidAlpha(2, 1) = 1.;
        break;
    case 3:
        dXidAlpha(2, 0) = 1.;
        dXidAlpha(1, 1) = 0.5;
        break;
    case 4:
        dXidAlpha(2, 0) = 1.;
        dXidAlpha(0, 1) = 0.5;
        dXidAlpha(1, 1) =-0.5;
        break;

    default:
        throw MechanicsException(__PRETTY_FUNCTION__, "PRISM3D has exactly five surfaces, 0 to 4. You tried to access " + std::to_string(rSurface) + ".");
    }
    return dXidAlpha;
}

std::vector<Eigen::VectorXd> NuTo::Interpolation3DPrism::GetSurfaceEdgesCoordinates(int rSurface) const
{
    std::vector<Eigen::VectorXd> surfaceEdgeCoordinates;
    switch (rSurface)
    {
    case 0:
        surfaceEdgeCoordinates.resize(3);
        surfaceEdgeCoordinates[0] = Eigen::Vector3d(0,0,-1);
        surfaceEdgeCoordinates[1] = Eigen::Vector3d(1,0,-1);
        surfaceEdgeCoordinates[2] = Eigen::Vector3d(0,1,-1);
        return surfaceEdgeCoordinates;
    case 1:
        surfaceEdgeCoordinates.resize(3);
        surfaceEdgeCoordinates[0] = Eigen::Vector3d(0,0, 1);
        surfaceEdgeCoordinates[1] = Eigen::Vector3d(1,0, 1);
        surfaceEdgeCoordinates[2] = Eigen::Vector3d(0,1, 1);
        return surfaceEdgeCoordinates;
    case 2:
        surfaceEdgeCoordinates.resize(4);
        surfaceEdgeCoordinates[0] = Eigen::Vector3d(0, 0,-1);
        surfaceEdgeCoordinates[1] = Eigen::Vector3d(1, 0,-1);
        surfaceEdgeCoordinates[2] = Eigen::Vector3d(1, 0, 1);
        surfaceEdgeCoordinates[3] = Eigen::Vector3d(0, 0, 1);
        return surfaceEdgeCoordinates;
    case 3:
        surfaceEdgeCoordinates.resize(4);
        surfaceEdgeCoordinates[0] = Eigen::Vector3d(0, 0,-1);
        surfaceEdgeCoordinates[1] = Eigen::Vector3d(0, 1,-1);
        surfaceEdgeCoordinates[2] = Eigen::Vector3d(0, 1, 1);
        surfaceEdgeCoordinates[3] = Eigen::Vector3d(0, 0, 1);
        return surfaceEdgeCoordinates;
    case 4:
        surfaceEdgeCoordinates.resize(4);
        surfaceEdgeCoordinates[0] = Eigen::Vector3d(0, 1,-1);
        surfaceEdgeCoordinates[1] = Eigen::Vector3d(1, 0,-1);
        surfaceEdgeCoordinates[2] = Eigen::Vector3d(1, 0, 1);
        surfaceEdgeCoordinates[3] = Eigen::Vector3d(0, 1, 1);
        return surfaceEdgeCoordinates;
    default:
        throw MechanicsException(__PRETTY_FUNCTION__, "PRISM3D has exactly five surfaces, 0 to 4. You tried to access " + std::to_string(rSurface) + ".");
    }
}
