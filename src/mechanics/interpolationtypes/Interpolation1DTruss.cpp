/*
 * Interpolation1DTruss.cpp
 *
 *  Created on: 8 May 2015
 *      Author: ttitsche
 */

#include "mechanics/MechanicsException.h"
#include "mechanics/elements/ElementShapeFunctions.h"
#include "mechanics/integrationtypes/IntegrationTypeEnum.h"
#include "mechanics/interpolationtypes/InterpolationTypeEnum.h"
#include "mechanics/interpolationtypes/Interpolation1DTruss.h"

NuTo::Interpolation1DTruss::Interpolation1DTruss(NuTo::Node::eDof rDofType, NuTo::Interpolation::eTypeOrder rTypeOrder,
                                                 int rDimension)
    : NuTo::Interpolation1D::Interpolation1D(rDofType, rTypeOrder, rDimension)
{
    Initialize();
}

NuTo::eIntegrationType NuTo::Interpolation1DTruss::GetStandardIntegrationType() const
{
    switch (mTypeOrder)
    {
    case NuTo::Interpolation::eTypeOrder::EQUIDISTANT1:
        return NuTo::eIntegrationType::IntegrationType1D2NGauss1Ip;
    case NuTo::Interpolation::eTypeOrder::EQUIDISTANT2:
        return NuTo::eIntegrationType::IntegrationType1D2NGauss2Ip;
    case NuTo::Interpolation::eTypeOrder::EQUIDISTANT3:
        return NuTo::eIntegrationType::IntegrationType1D2NGauss3Ip;
    case NuTo::Interpolation::eTypeOrder::EQUIDISTANT4:
        return NuTo::eIntegrationType::IntegrationType1D2NGauss4Ip;
    case NuTo::Interpolation::eTypeOrder::LOBATTO2:
        return NuTo::eIntegrationType::IntegrationType1D2NLobatto3Ip;
    case NuTo::Interpolation::eTypeOrder::LOBATTO3:
        return NuTo::eIntegrationType::IntegrationType1D2NLobatto4Ip;
    case NuTo::Interpolation::eTypeOrder::LOBATTO4:
        return NuTo::eIntegrationType::IntegrationType1D2NLobatto5Ip;
    default:
        throw MechanicsException(__PRETTY_FUNCTION__, "Interpolation for exact integration of " +
                                                              Interpolation::TypeOrderToString(mTypeOrder) +
                                                              " not implemented");
    }
}

Eigen::VectorXd NuTo::Interpolation1DTruss::CalculateNaturalNodeCoordinates(int rNodeIndexDof) const
{
    switch (mTypeOrder)
    {
    case NuTo::Interpolation::eTypeOrder::EQUIDISTANT1:
        return ShapeFunctions1D::NodeCoordinatesTrussOrder1(rNodeIndexDof);
    case NuTo::Interpolation::eTypeOrder::EQUIDISTANT2:
    case NuTo::Interpolation::eTypeOrder::LOBATTO2:
        return ShapeFunctions1D::NodeCoordinatesTrussOrder2(rNodeIndexDof);
    case NuTo::Interpolation::eTypeOrder::EQUIDISTANT3:
        return ShapeFunctions1D::NodeCoordinatesTrussOrder3(rNodeIndexDof);
    case NuTo::Interpolation::eTypeOrder::EQUIDISTANT4:
        return ShapeFunctions1D::NodeCoordinatesTrussOrder4(rNodeIndexDof);
    case NuTo::Interpolation::eTypeOrder::LOBATTO3:
        return ShapeFunctions1D::NodeCoordinatesTrussSpectralOrder3(rNodeIndexDof);
    case NuTo::Interpolation::eTypeOrder::LOBATTO4:
        return ShapeFunctions1D::NodeCoordinatesTrussSpectralOrder4(rNodeIndexDof);
    default:
        throw MechanicsException(__PRETTY_FUNCTION__, "Node arrangement for " +
                                                              Interpolation::TypeOrderToString(mTypeOrder) +
                                                              " not implemented");
    }
}

Eigen::VectorXd NuTo::Interpolation1DTruss::CalculateShapeFunctions(const Eigen::VectorXd& rCoordinates) const
{
    switch (mTypeOrder)
    {
    case NuTo::Interpolation::eTypeOrder::EQUIDISTANT1:
        return ShapeFunctions1D::ShapeFunctionsTrussOrder1(rCoordinates);
    case NuTo::Interpolation::eTypeOrder::EQUIDISTANT2:
    case NuTo::Interpolation::eTypeOrder::LOBATTO2:
        return ShapeFunctions1D::ShapeFunctionsTrussOrder2(rCoordinates);
    case NuTo::Interpolation::eTypeOrder::EQUIDISTANT3:
        return ShapeFunctions1D::ShapeFunctionsTrussOrder3(rCoordinates);
    case NuTo::Interpolation::eTypeOrder::EQUIDISTANT4:
        return ShapeFunctions1D::ShapeFunctionsTrussOrder4(rCoordinates);
    case NuTo::Interpolation::eTypeOrder::LOBATTO3:
        return ShapeFunctions1D::ShapeFunctionsTrussSpectralOrder3(rCoordinates);
    case NuTo::Interpolation::eTypeOrder::LOBATTO4:
        return ShapeFunctions1D::ShapeFunctionsTrussSpectralOrder4(rCoordinates);
    default:
        throw MechanicsException(__PRETTY_FUNCTION__, "Interpolation order for " +
                                                              Interpolation::TypeOrderToString(mTypeOrder) +
                                                              " not implemented");
    }
}

Eigen::MatrixXd
NuTo::Interpolation1DTruss::CalculateDerivativeShapeFunctionsNatural(const Eigen::VectorXd& rCoordinates) const
{
    switch (mTypeOrder)
    {
    case NuTo::Interpolation::eTypeOrder::EQUIDISTANT1:
        return ShapeFunctions1D::DerivativeShapeFunctionsTrussOrder1();
    case NuTo::Interpolation::eTypeOrder::EQUIDISTANT2:
    case NuTo::Interpolation::eTypeOrder::LOBATTO2:
        return ShapeFunctions1D::DerivativeShapeFunctionsTrussOrder2(rCoordinates);
    case NuTo::Interpolation::eTypeOrder::EQUIDISTANT3:
        return ShapeFunctions1D::DerivativeShapeFunctionsTrussOrder3(rCoordinates);
    case NuTo::Interpolation::eTypeOrder::EQUIDISTANT4:
        return ShapeFunctions1D::DerivativeShapeFunctionsTrussOrder4(rCoordinates);
    case NuTo::Interpolation::eTypeOrder::LOBATTO3:
        return ShapeFunctions1D::DerivativeShapeFunctionsTrussSpectralOrder3(rCoordinates);
    case NuTo::Interpolation::eTypeOrder::LOBATTO4:
        return ShapeFunctions1D::DerivativeShapeFunctionsTrussSpectralOrder4(rCoordinates);
    default:
        throw MechanicsException(__PRETTY_FUNCTION__, "Interpolation order for " +
                                                              Interpolation::TypeOrderToString(mTypeOrder) +
                                                              " not implemented");
    }
}

Eigen::VectorXd NuTo::Interpolation1DTruss::CalculateNaturalSurfaceCoordinates(const Eigen::VectorXd&,
                                                                               int rSurface) const
{
    Eigen::VectorXd naturalCoordinates(1);
    switch (rSurface)
    {
    case 0:
        naturalCoordinates(0, 0) = -1;
        break;
    case 1:
        naturalCoordinates(0, 0) = 1;
        break;
    default:
        throw MechanicsException(__PRETTY_FUNCTION__, "TRUSS1D has exactly two surfaces, 0 or 1. You tried to access " +
                                                              std::to_string(rSurface) + ".");
    }
    return naturalCoordinates;
}

Eigen::MatrixXd NuTo::Interpolation1DTruss::CalculateDerivativeNaturalSurfaceCoordinates(const Eigen::VectorXd&,
                                                                                         int) const
{
    Eigen::VectorXd naturalSurfaceCoordinates(1);
    // calculating a transformation from 1D --> 0D returns 0
    // Since we are in 3D, the boundary turns into a surface with the elements area
    // Thus, returning 1 here.
    naturalSurfaceCoordinates(0, 0) = 1;
    return naturalSurfaceCoordinates;
}

int NuTo::Interpolation1DTruss::CalculateNumNodes() const
{
    switch (mTypeOrder)
    {
    case NuTo::Interpolation::eTypeOrder::EQUIDISTANT1:
        return 2;
    case NuTo::Interpolation::eTypeOrder::EQUIDISTANT2:
    case NuTo::Interpolation::eTypeOrder::LOBATTO2:
        return 3;
    case NuTo::Interpolation::eTypeOrder::EQUIDISTANT3:
    case NuTo::Interpolation::eTypeOrder::LOBATTO3:
        return 4;
    case NuTo::Interpolation::eTypeOrder::EQUIDISTANT4:
    case NuTo::Interpolation::eTypeOrder::LOBATTO4:
        return 5;
    default:
        throw MechanicsException(__PRETTY_FUNCTION__, "Interpolation type and order " +
                                                              Interpolation::TypeOrderToString(mTypeOrder) +
                                                              " not implemented");
    }
}
