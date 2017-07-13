/*
 * InterpolationType2DQuad.cpp
 *
 *  Created on: 24 Apr 2015
 *      Author: ttitsche
 */

#include "base/Exception.h"
#include "mechanics/elements/ElementShapeFunctions.h"
#include "mechanics/integrationtypes/IntegrationTypeEnum.h"
#include "mechanics/interpolationtypes/InterpolationTypeEnum.h"
#include "mechanics/interpolationtypes/Interpolation2DQuad.h"

NuTo::Interpolation2DQuad::Interpolation2DQuad(NuTo::Node::eDof rDofType, NuTo::Interpolation::eTypeOrder rTypeOrder,
                                               int rDimension)
    : Interpolation2D::Interpolation2D(rDofType, rTypeOrder, rDimension)
{
    Initialize();
}

NuTo::eIntegrationType NuTo::Interpolation2DQuad::GetStandardIntegrationType() const
{
    switch (mTypeOrder)
    {
    case NuTo::Interpolation::eTypeOrder::EQUIDISTANT1:
        return NuTo::eIntegrationType::IntegrationType2D4NGauss4Ip;
    case NuTo::Interpolation::eTypeOrder::EQUIDISTANT2:
        return NuTo::eIntegrationType::IntegrationType2D4NGauss4Ip;
    case NuTo::Interpolation::eTypeOrder::LOBATTO2:
        return NuTo::eIntegrationType::IntegrationType2D4NLobatto9Ip;
    case NuTo::Interpolation::eTypeOrder::LOBATTO3:
        return NuTo::eIntegrationType::IntegrationType2D4NLobatto16Ip;
    case NuTo::Interpolation::eTypeOrder::LOBATTO4:
        return NuTo::eIntegrationType::IntegrationType2D4NLobatto25Ip;

    default:
        throw Exception(__PRETTY_FUNCTION__, "Interpolation for exact integration of " +
                                                              Interpolation::TypeOrderToString(mTypeOrder) +
                                                              " not implemented");
    }
}

Eigen::VectorXd NuTo::Interpolation2DQuad::CalculateShapeFunctions(const Eigen::VectorXd& rCoordinates) const
{
    switch (mTypeOrder)
    {
    case NuTo::Interpolation::eTypeOrder::EQUIDISTANT1:
        return ShapeFunctions2D::ShapeFunctionsQuadOrder1(rCoordinates);
    case NuTo::Interpolation::eTypeOrder::EQUIDISTANT2:
        return ShapeFunctions2D::ShapeFunctionsQuadOrder2(rCoordinates);
    case NuTo::Interpolation::eTypeOrder::LOBATTO2:
        return ShapeFunctions2D::ShapeFunctionsQuadSpectralOrder2(rCoordinates);
    case NuTo::Interpolation::eTypeOrder::LOBATTO3:
        return ShapeFunctions2D::ShapeFunctionsQuadSpectralOrder3(rCoordinates);
    case NuTo::Interpolation::eTypeOrder::LOBATTO4:
        return ShapeFunctions2D::ShapeFunctionsQuadSpectralOrder4(rCoordinates);
    default:
        throw Exception(__PRETTY_FUNCTION__, "Interpolation order for " +
                                                              Interpolation::TypeOrderToString(mTypeOrder) +
                                                              " not implemented");
    }
}

Eigen::MatrixXd
NuTo::Interpolation2DQuad::CalculateDerivativeShapeFunctionsNatural(const Eigen::VectorXd& rCoordinates) const
{
    switch (mTypeOrder)
    {
    case NuTo::Interpolation::eTypeOrder::EQUIDISTANT1:
        return ShapeFunctions2D::DerivativeShapeFunctionsQuadOrder1(rCoordinates);
    case NuTo::Interpolation::eTypeOrder::EQUIDISTANT2:
        return ShapeFunctions2D::DerivativeShapeFunctionsQuadOrder2(rCoordinates);
    case NuTo::Interpolation::eTypeOrder::LOBATTO2:
        return ShapeFunctions2D::DerivativeShapeFunctionsQuadSpectralOrder2(rCoordinates);
    case NuTo::Interpolation::eTypeOrder::LOBATTO3:
        return ShapeFunctions2D::DerivativeShapeFunctionsQuadSpectralOrder3(rCoordinates);
    case NuTo::Interpolation::eTypeOrder::LOBATTO4:
        return ShapeFunctions2D::DerivativeShapeFunctionsQuadSpectralOrder4(rCoordinates);

    default:
        throw Exception(__PRETTY_FUNCTION__, "Interpolation order for " +
                                                              Interpolation::TypeOrderToString(mTypeOrder) +
                                                              " not implemented");
    }
}

Eigen::VectorXd NuTo::Interpolation2DQuad::CalculateNaturalNodeCoordinates(int rNodeIndexDof) const
{
    switch (mTypeOrder)
    {
    case NuTo::Interpolation::eTypeOrder::EQUIDISTANT1:
        return ShapeFunctions2D::NodeCoordinatesQuadOrder1(rNodeIndexDof);
    case NuTo::Interpolation::eTypeOrder::EQUIDISTANT2:
        return ShapeFunctions2D::NodeCoordinatesQuadOrder2(rNodeIndexDof);
    case NuTo::Interpolation::eTypeOrder::LOBATTO2:
        return ShapeFunctions2D::NodeCoordinatesQuadSpectralOrder2(rNodeIndexDof);
    case NuTo::Interpolation::eTypeOrder::LOBATTO3:
        return ShapeFunctions2D::NodeCoordinatesQuadSpectralOrder3(rNodeIndexDof);
    case NuTo::Interpolation::eTypeOrder::LOBATTO4:
        return ShapeFunctions2D::NodeCoordinatesQuadSpectralOrder4(rNodeIndexDof);

    default:
        throw Exception(__PRETTY_FUNCTION__, "Node arrangement for " +
                                                              Interpolation::TypeOrderToString(mTypeOrder) +
                                                              " not implemented");
    }
}

Eigen::VectorXd
NuTo::Interpolation2DQuad::CalculateNaturalSurfaceCoordinates(const Eigen::VectorXd& rNaturalSurfaceCoordinates,
                                                              int rSurface) const
{
    assert(rNaturalSurfaceCoordinates.rows() == 1);
    switch (rSurface)
    {
    case 0:
        return Eigen::Vector2d(rNaturalSurfaceCoordinates(0), -1.);
    case 1:
        return Eigen::Vector2d(1., rNaturalSurfaceCoordinates(0));
    case 2:
        return Eigen::Vector2d(-rNaturalSurfaceCoordinates(0), 1.);
    case 3:
        return Eigen::Vector2d(-1., -rNaturalSurfaceCoordinates(0));
    default:
        throw Exception(__PRETTY_FUNCTION__, "QUAD2D has exactly four surfaces, 0 to 3. You tried to access " +
                                                              std::to_string(rSurface) + ".");
    }
}

Eigen::MatrixXd NuTo::Interpolation2DQuad::CalculateDerivativeNaturalSurfaceCoordinates(
        const Eigen::VectorXd& rNaturalSurfaceCoordinates, int rSurface) const
{
    assert(rNaturalSurfaceCoordinates.rows() == 1);
    switch (rSurface)
    {
    case 0:
        return Eigen::Vector2d(1, 0);
    case 1:
        return Eigen::Vector2d(0, 1);
    case 2:
        return Eigen::Vector2d(-1, 0);
    case 3:
        return Eigen::Vector2d(0, -1);
    default:
        throw Exception(__PRETTY_FUNCTION__, "QUAD2D has exactly four surfaces, 0 to 3. You tried to access " +
                                                              std::to_string(rSurface) + ".");
    }
}

int NuTo::Interpolation2DQuad::CalculateNumNodes() const
{
    switch (mTypeOrder)
    {
    case NuTo::Interpolation::eTypeOrder::EQUIDISTANT1:
        return 4;
    case NuTo::Interpolation::eTypeOrder::EQUIDISTANT2:
        return 8;
    case NuTo::Interpolation::eTypeOrder::LOBATTO2:
        return 9;
    case NuTo::Interpolation::eTypeOrder::LOBATTO3:
        return 16;
    case NuTo::Interpolation::eTypeOrder::LOBATTO4:
        return 25;

    default:
        throw Exception(__PRETTY_FUNCTION__, "Interpolation type and order " +
                                                              Interpolation::TypeOrderToString(mTypeOrder) +
                                                              " not implemented");
    }
}
