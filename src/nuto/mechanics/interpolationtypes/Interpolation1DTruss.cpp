/*
 * Interpolation1DTruss.cpp
 *
 *  Created on: 8 May 2015
 *      Author: ttitsche
 */

#include "nuto/mechanics/interpolationtypes/Interpolation1DTruss.h"

NuTo::Interpolation1DTruss::Interpolation1DTruss(NuTo::Node::eAttributes rDofType, NuTo::Interpolation::eTypeOrder rTypeOrder, int rDimension) :
        NuTo::Interpolation1D::Interpolation1D(rDofType, rTypeOrder, rDimension)
{
    Initialize();
}

NuTo::IntegrationType::eIntegrationType NuTo::Interpolation1DTruss::GetStandardIntegrationType() const
{
    switch (mTypeOrder)
    {
    case NuTo::Interpolation::eTypeOrder::EQUIDISTANT1:
        return NuTo::IntegrationType::IntegrationType1D2NGauss1Ip;
    case NuTo::Interpolation::eTypeOrder::EQUIDISTANT2:
        return NuTo::IntegrationType::IntegrationType1D2NGauss2Ip;
    case NuTo::Interpolation::eTypeOrder::EQUIDISTANT3:
        return NuTo::IntegrationType::IntegrationType1D2NGauss3Ip;
    case NuTo::Interpolation::eTypeOrder::EQUIDISTANT4:
            return NuTo::IntegrationType::IntegrationType1D2NGauss4Ip;
    case NuTo::Interpolation::eTypeOrder::LOBATTO2:
        return NuTo::IntegrationType::IntegrationType1D2NLobatto3Ip;
    case NuTo::Interpolation::eTypeOrder::LOBATTO3:
        return NuTo::IntegrationType::IntegrationType1D2NLobatto4Ip;
    case NuTo::Interpolation::eTypeOrder::LOBATTO4:
        return NuTo::IntegrationType::IntegrationType1D2NLobatto5Ip;
    default:
        throw MechanicsException("[NuTo::Interpolation1DTruss::GetStandardIntegrationType] Interpolation for exact integration of " + Interpolation::TypeOrderToString(mTypeOrder) + " not implemented");
    }
    return NuTo::IntegrationType::IntegrationType1D2NGauss2Ip;
}

const Eigen::VectorXd NuTo::Interpolation1DTruss::CalculateNaturalNodeCoordinates(int rNodeIndexDof) const
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
        throw MechanicsException("[NuTo::Interpolation1DTruss::CalculateNaturalNodeCoordinates] Node arrangement for " + Interpolation::TypeOrderToString(mTypeOrder) + " not implemented");
    }
}

const Eigen::VectorXd NuTo::Interpolation1DTruss::CalculateShapeFunctions(const Eigen::VectorXd& rCoordinates) const
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
        throw MechanicsException("[NuTo::Interpolation1DTruss::CalculateShapeFunctions] Interpolation order for " + Interpolation::TypeOrderToString(mTypeOrder) + " not implemented");
    }
}

const Eigen::MatrixXd NuTo::Interpolation1DTruss::CalculateDerivativeShapeFunctionsNatural(const Eigen::VectorXd& rCoordinates) const
{
    switch (mTypeOrder)
    {
    case NuTo::Interpolation::eTypeOrder::EQUIDISTANT1:
        return ShapeFunctions1D::DerivativeShapeFunctionsTrussOrder1(rCoordinates);
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
        throw MechanicsException("[NuTo::Interpolation1DTruss::CalculateDerivativeShapeFunctionsNatural] Interpolation order for " + Interpolation::TypeOrderToString(mTypeOrder) + " not implemented");
    }
}

const Eigen::VectorXd NuTo::Interpolation1DTruss::CalculateNaturalSurfaceCoordinates(const Eigen::VectorXd& rNaturalSurfaceCoordinates, int rSurface) const
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
        throw("[NuTo::Interpolation1DTruss::CalculateNaturalSurfaceCoordinates] TRUSS1D has exactly two surfaces, 0 or 1. You tried to access " + std::to_string(rSurface) + ".");
    }
    return naturalCoordinates;
}

const Eigen::MatrixXd NuTo::Interpolation1DTruss::CalculateDerivativeNaturalSurfaceCoordinates(const Eigen::VectorXd& rNaturalSurfaceCoordinates, int rSurface) const
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
        throw MechanicsException("[NuTo::Interpolation1DTruss::GetNumNodes] Interpolation type and order " + Interpolation::TypeOrderToString(mTypeOrder) + " not implemented");
    }

}

#ifdef ENABLE_SERIALIZATION
template<class Archive>
void NuTo::Interpolation1DTruss::serialize(Archive & ar, const unsigned int version)
{
#ifdef DEBUG_SERIALIZATION
    std::cout << "start serialize Interpolation1D" << std::endl;
#endif
    ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(Interpolation1D);
#ifdef DEBUG_SERIALIZATION
    std::cout << "finish serialize Interpolation1D" << std::endl;
#endif
}
BOOST_CLASS_EXPORT_IMPLEMENT(NuTo::Interpolation1DTruss)
#endif  // ENABLE_SERIALIZATION
