
#include "mechanics/MechanicsException.h"
#include "mechanics/elements/ElementShapeFunctions.h"
#include "mechanics/integrationtypes/IntegrationTypeEnum.h"
#include "mechanics/interpolationtypes/InterpolationTypeEnum.h"
#include "mechanics/interpolationtypes/Interpolation1DInterface.h"

NuTo::Interpolation1DInterface::Interpolation1DInterface(NuTo::Node::eDof rDofType, NuTo::Interpolation::eTypeOrder rTypeOrder, int rDimension) :
        NuTo::Interpolation1DTruss::Interpolation1DTruss(rDofType, rTypeOrder, rDimension)

{
    Initialize();
}

NuTo::eIntegrationType NuTo::Interpolation1DInterface::GetStandardIntegrationType() const
{
    switch (mTypeOrder)
    {
    case NuTo::Interpolation::eTypeOrder::EQUIDISTANT1:
    {
//        return NuTo::eIntegrationType::IntegrationType1D2NGauss2Ip;
        return NuTo::eIntegrationType::IntegrationType1D2NLobatto3Ip;
    }
    case NuTo::Interpolation::eTypeOrder::EQUIDISTANT2:
    {
        return NuTo::eIntegrationType::IntegrationType1D2NGauss3Ip;
    }
    default:
        throw MechanicsException(__PRETTY_FUNCTION__, "Interpolation for exact integration of " + Interpolation::TypeOrderToString(mTypeOrder) + " not implemented");
    }
}

Eigen::VectorXd NuTo::Interpolation1DInterface::CalculateNaturalNodeCoordinates(int rNodeIndexDof) const
{
    switch (mTypeOrder)
    {
    case NuTo::Interpolation::eTypeOrder::EQUIDISTANT1:
        return ShapeFunctionsInterface2D::NodeCoordinatesInterface2dOrder1(rNodeIndexDof);
    case NuTo::Interpolation::eTypeOrder::EQUIDISTANT2:
        return ShapeFunctionsInterface2D::NodeCoordinatesInterface2dOrder2(rNodeIndexDof);
    default:
        throw MechanicsException(__PRETTY_FUNCTION__, "Node arrangement for " + Interpolation::TypeOrderToString(mTypeOrder) + " not implemented");
    }
}

Eigen::VectorXd NuTo::Interpolation1DInterface::CalculateShapeFunctions(const Eigen::VectorXd& rCoordinates) const
{
    switch (mTypeOrder)
    {
    case NuTo::Interpolation::eTypeOrder::EQUIDISTANT1:
        return ShapeFunctionsInterface2D::ShapeFunctionsInterface2dOrder1(rCoordinates);
    case NuTo::Interpolation::eTypeOrder::EQUIDISTANT2:
        return ShapeFunctionsInterface2D::ShapeFunctionsInterface2dOrder2(rCoordinates);
    default:
        throw MechanicsException(__PRETTY_FUNCTION__, "Interpolation order for " + Interpolation::TypeOrderToString(mTypeOrder) + " not implemented");
    }
}

Eigen::MatrixXd NuTo::Interpolation1DInterface::CalculateDerivativeShapeFunctionsNatural(const Eigen::VectorXd& rCoordinates) const
{
    switch (mTypeOrder)
    {
    case NuTo::Interpolation::eTypeOrder::EQUIDISTANT1:
        return ShapeFunctionsInterface2D::DerivativeShapeFunctionsInterface2dOrder1(rCoordinates);
    case NuTo::Interpolation::eTypeOrder::EQUIDISTANT2:
        return ShapeFunctionsInterface2D::DerivativeShapeFunctionsInterface2dOrder2(rCoordinates);
    default:
        throw MechanicsException(__PRETTY_FUNCTION__, "Interpolation order for " + Interpolation::TypeOrderToString(mTypeOrder) + " not implemented");
    }
}

Eigen::VectorXd NuTo::Interpolation1DInterface::CalculateNaturalSurfaceCoordinates(const Eigen::VectorXd& rNaturalSurfaceCoordinates, int rSurface) const
{
    throw("[NuTo::Interpolation1DInterface::CalculateNaturalSurfaceCoordinates] not implemented");
}


Eigen::MatrixXd NuTo::Interpolation1DInterface::CalculateDerivativeNaturalSurfaceCoordinates(const Eigen::VectorXd& rNaturalSurfaceCoordinates, int rSurface) const
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
        throw MechanicsException(__PRETTY_FUNCTION__, "Interpolation type and order " + Interpolation::TypeOrderToString(mTypeOrder) + " not implemented");
    }

}


#ifdef ENABLE_SERIALIZATION
BOOST_CLASS_EXPORT_IMPLEMENT(NuTo::Interpolation1DInterface)
#endif

