
#include "nuto/mechanics/interpolationtypes/Interpolation1DInterface.h"

NuTo::Interpolation1DInterface::Interpolation1DInterface(NuTo::Node::eAttributes rDofType, NuTo::Interpolation::eTypeOrder rTypeOrder, int rDimension) :
        NuTo::Interpolation1DTruss::Interpolation1DTruss(rDofType, rTypeOrder, rDimension)
{
    Initialize();
}

NuTo::IntegrationType::eIntegrationType NuTo::Interpolation1DInterface::GetStandardIntegrationType() const
{
    switch (mTypeOrder)
    {
    case NuTo::Interpolation::eTypeOrder::EQUIDISTANT1:
    {
//        return NuTo::IntegrationType::IntegrationType1D2NGauss2Ip;
        return NuTo::IntegrationType::IntegrationType1D2NLobatto3Ip;
    }
    case NuTo::Interpolation::eTypeOrder::EQUIDISTANT2:
    {
        return NuTo::IntegrationType::IntegrationType1D2NGauss3Ip;
    }
    default:
        throw MechanicsException("[NuTo::Interpolation1DInterface::GetStandardIntegrationType] Interpolation for exact integration of " + Interpolation::TypeOrderToString(mTypeOrder) + " not implemented");
    }
}

const Eigen::VectorXd NuTo::Interpolation1DInterface::CalculateNaturalNodeCoordinates(int rNodeIndexDof) const
{
    switch (mTypeOrder)
    {
    case NuTo::Interpolation::eTypeOrder::EQUIDISTANT1:
        return ShapeFunctionsInterface2D::NodeCoordinatesInterface2dOrder1(rNodeIndexDof);
    case NuTo::Interpolation::eTypeOrder::EQUIDISTANT2:
        return ShapeFunctionsInterface2D::NodeCoordinatesInterface2dOrder2(rNodeIndexDof);
    default:
        throw MechanicsException("[NuTo::Interpolation1DInterface::CalculateNaturalNodeCoordinates] Node arrangement for " + Interpolation::TypeOrderToString(mTypeOrder) + " not implemented");
    }
}

const Eigen::VectorXd NuTo::Interpolation1DInterface::CalculateShapeFunctions(const Eigen::VectorXd& rCoordinates) const
{
    switch (mTypeOrder)
    {
    case NuTo::Interpolation::eTypeOrder::EQUIDISTANT1:
        return ShapeFunctionsInterface2D::ShapeFunctionsInterface2dOrder1(rCoordinates);
    case NuTo::Interpolation::eTypeOrder::EQUIDISTANT2:
        return ShapeFunctionsInterface2D::ShapeFunctionsInterface2dOrder2(rCoordinates);
    default:
        throw MechanicsException("[NuTo::Interpolation1DInterface::CalculateShapeFunctions] Interpolation order for " + Interpolation::TypeOrderToString(mTypeOrder) + " not implemented");
    }
}

const Eigen::MatrixXd NuTo::Interpolation1DInterface::CalculateDerivativeShapeFunctionsNatural(const Eigen::VectorXd& rCoordinates) const
{
    switch (mTypeOrder)
    {
    case NuTo::Interpolation::eTypeOrder::EQUIDISTANT1:
        return ShapeFunctionsInterface2D::DerivativeShapeFunctionsInterface2dOrder1(rCoordinates);
    case NuTo::Interpolation::eTypeOrder::EQUIDISTANT2:
        return ShapeFunctionsInterface2D::DerivativeShapeFunctionsInterface2dOrder2(rCoordinates);
    default:
        throw MechanicsException("[NuTo::Interpolation1DInterface::CalculateDerivativeShapeFunctionsNatural] Interpolation order for " + Interpolation::TypeOrderToString(mTypeOrder) + " not implemented");
    }
}

const Eigen::VectorXd NuTo::Interpolation1DInterface::CalculateNaturalSurfaceCoordinates(const Eigen::VectorXd& rNaturalSurfaceCoordinates, int rSurface) const
{
    throw("[NuTo::Interpolation1DInterface::CalculateNaturalSurfaceCoordinates] not implemented");
}

const Eigen::MatrixXd NuTo::Interpolation1DInterface::CalculateDerivativeNaturalSurfaceCoordinates(const Eigen::VectorXd& rNaturalSurfaceCoordinates, int rSurface) const
{
    throw("[NuTo::Interpolation1DInterface::CalculateDerivativeNaturalSurfaceCoordinates] not implemented");
}

int NuTo::Interpolation1DInterface::CalculateNumNodes() const
{
    switch (mTypeOrder)
    {
    case NuTo::Interpolation::eTypeOrder::EQUIDISTANT1:
        return 4;
    case NuTo::Interpolation::eTypeOrder::EQUIDISTANT2:
        return 6;
    default:
        throw MechanicsException("[NuTo::Interpolation1DInterface::GetNumNodes] Interpolation type and order " + Interpolation::TypeOrderToString(mTypeOrder) + " not implemented");
    }

}
