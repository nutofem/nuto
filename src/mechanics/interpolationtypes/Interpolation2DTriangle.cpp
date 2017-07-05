/*
 * Interpolation2DTriangle.cpp
 *
 *  Created on: 20 Mar 2015
 *      Author: ttitsche
 */

#include "base/Exception.h"
#include "mechanics/elements/ElementShapeFunctions.h"
#include "mechanics/integrationtypes/IntegrationTypeEnum.h"
#include "mechanics/interpolationtypes/InterpolationTypeEnum.h"
#include "mechanics/interpolationtypes/Interpolation2DTriangle.h"

NuTo::Interpolation2DTriangle::Interpolation2DTriangle(NuTo::Node::eDof rDofType,
                                                       NuTo::Interpolation::eTypeOrder rTypeOrder, int rDimension)
    : Interpolation2D::Interpolation2D(rDofType, rTypeOrder, rDimension)
{
    Initialize();
}

NuTo::eIntegrationType NuTo::Interpolation2DTriangle::GetStandardIntegrationType() const
{
    switch (mTypeOrder)
    {
    case NuTo::Interpolation::eTypeOrder::EQUIDISTANT1:
        return NuTo::eIntegrationType::IntegrationType2D3NGauss1Ip;
    case NuTo::Interpolation::eTypeOrder::EQUIDISTANT2:
        return NuTo::eIntegrationType::IntegrationType2D3NGauss3Ip;
    case NuTo::Interpolation::eTypeOrder::EQUIDISTANT3:
        return NuTo::eIntegrationType::IntegrationType2D3NGauss6Ip;
    case NuTo::Interpolation::eTypeOrder::EQUIDISTANT4:
        return NuTo::eIntegrationType::IntegrationType2D3NGauss12Ip;
    default:
        throw Exception(__PRETTY_FUNCTION__, "Interpolation for exact integration of " + Interpolation::TypeOrderToString(mTypeOrder) + " not implemented");
    }
}

Eigen::VectorXd NuTo::Interpolation2DTriangle::CalculateShapeFunctions(const Eigen::VectorXd& rCoordinates) const
{
    switch (mTypeOrder)
    {
    case NuTo::Interpolation::eTypeOrder::EQUIDISTANT1:
        return ShapeFunctions2D::ShapeFunctionsTriangleOrder1(rCoordinates);
    case NuTo::Interpolation::eTypeOrder::EQUIDISTANT2:
        return ShapeFunctions2D::ShapeFunctionsTriangleOrder2(rCoordinates);
    case NuTo::Interpolation::eTypeOrder::EQUIDISTANT3:
        return ShapeFunctions2D::ShapeFunctionsTriangleOrder3(rCoordinates);
    case NuTo::Interpolation::eTypeOrder::EQUIDISTANT4:
        return ShapeFunctions2D::ShapeFunctionsTriangleOrder4(rCoordinates);
    default:
        throw Exception(__PRETTY_FUNCTION__, "Interpolation order for " + Interpolation::TypeOrderToString(mTypeOrder) + " not implemented");
    }
}

Eigen::MatrixXd
NuTo::Interpolation2DTriangle::CalculateDerivativeShapeFunctionsNatural(const Eigen::VectorXd& rCoordinates) const
{
    switch (mTypeOrder)
    {
    case NuTo::Interpolation::eTypeOrder::EQUIDISTANT1:
        return ShapeFunctions2D::DerivativeShapeFunctionsTriangleOrder1(rCoordinates);
    case NuTo::Interpolation::eTypeOrder::EQUIDISTANT2:
        return ShapeFunctions2D::DerivativeShapeFunctionsTriangleOrder2(rCoordinates);
    case NuTo::Interpolation::eTypeOrder::EQUIDISTANT3:
        return ShapeFunctions2D::DerivativeShapeFunctionsTriangleOrder3(rCoordinates);
    case NuTo::Interpolation::eTypeOrder::EQUIDISTANT4:
        return ShapeFunctions2D::DerivativeShapeFunctionsTriangleOrder4(rCoordinates);
    default:
        throw Exception(__PRETTY_FUNCTION__, "Interpolation order for " + Interpolation::TypeOrderToString(mTypeOrder) + " not implemented");
    }
}

Eigen::VectorXd NuTo::Interpolation2DTriangle::CalculateNaturalNodeCoordinates(int rNodeIndexDof) const
{
    switch (mTypeOrder)
    {
    case NuTo::Interpolation::eTypeOrder::EQUIDISTANT1:
        return ShapeFunctions2D::NodeCoordinatesTriangleOrder1(rNodeIndexDof);
    case NuTo::Interpolation::eTypeOrder::EQUIDISTANT2:
        return ShapeFunctions2D::NodeCoordinatesTriangleOrder2(rNodeIndexDof);
    case NuTo::Interpolation::eTypeOrder::EQUIDISTANT3:
        return ShapeFunctions2D::NodeCoordinatesTriangleOrder3(rNodeIndexDof);
    case NuTo::Interpolation::eTypeOrder::EQUIDISTANT4:
        return ShapeFunctions2D::NodeCoordinatesTriangleOrder4(rNodeIndexDof);
    default:
        throw Exception(__PRETTY_FUNCTION__, "Node arrangement for " + Interpolation::TypeOrderToString(mTypeOrder) + " not implemented");
    }
}

int NuTo::Interpolation2DTriangle::CalculateNumNodes() const
{
    // NumNodes(i) = SUM_0^i (i+1)
    switch (mTypeOrder)
    {
    case NuTo::Interpolation::eTypeOrder::EQUIDISTANT1:
        return 3;
    case NuTo::Interpolation::eTypeOrder::EQUIDISTANT2:
        return 6;
    case NuTo::Interpolation::eTypeOrder::EQUIDISTANT3:
        return 10;
    case NuTo::Interpolation::eTypeOrder::EQUIDISTANT4:
        return 15;
    default:
        throw Exception(__PRETTY_FUNCTION__, "Interpolation type and order " + Interpolation::TypeOrderToString(mTypeOrder) + " not implemented");
    }
}

Eigen::VectorXd
NuTo::Interpolation2DTriangle::CalculateNaturalSurfaceCoordinates(const Eigen::VectorXd& rNaturalSurfaceCoordinates,
                                                                  int rSurface) const
{
    assert(rNaturalSurfaceCoordinates.rows() == 1);
    switch (rSurface)
    {
    case 0:
        return Eigen::Vector2d(.5 * (1 + rNaturalSurfaceCoordinates(0)), 0);
    case 1:
        return Eigen::Vector2d(.5 * (1 - rNaturalSurfaceCoordinates(0)), .5 * (1 + rNaturalSurfaceCoordinates(0)));
    case 2:
        return Eigen::Vector2d(0, .5 * (1 - rNaturalSurfaceCoordinates(0)));
    default:
        throw Exception(__PRETTY_FUNCTION__, "TRIANGLE2D has exactly three surfaces, 0 to 2. You tried to access " + std::to_string(rSurface) + ".");
    }
}

Eigen::MatrixXd NuTo::Interpolation2DTriangle::CalculateDerivativeNaturalSurfaceCoordinates(
        const Eigen::VectorXd& rNaturalSurfaceCoordinates, int rSurface) const
{
    assert(rNaturalSurfaceCoordinates.rows() == 1);
    switch (rSurface)
    {
    case 0:
        return Eigen::Vector2d(.5, 0);
    case 1:
        return Eigen::Vector2d(-.5, .5);
    case 2:
        return Eigen::Vector2d(0, -.5);
    default:
        throw Exception(__PRETTY_FUNCTION__, "TRIANGLE2D has exactly three surfaces, 0 to 2. You tried to access " + std::to_string(rSurface) + ".");
    }
}
