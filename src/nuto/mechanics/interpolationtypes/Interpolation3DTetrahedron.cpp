/*
 * Interpolation3DTetrahedron.h
 *
 *  Created on: 19 May 2015
 *      Author: ttitsche
 */

#include "nuto/mechanics/MechanicsException.h"
#include "nuto/mechanics/elements/ElementShapeFunctions.h"
#include "nuto/mechanics/integrationtypes/IntegrationTypeEnum.h"
#include "nuto/mechanics/interpolationtypes/InterpolationTypeEnum.h"
#include "nuto/mechanics/interpolationtypes/Interpolation3DTetrahedron.h"

NuTo::Interpolation3DTetrahedron::Interpolation3DTetrahedron(NuTo::Node::eDof rDofType, NuTo::Interpolation::eTypeOrder rTypeOrder, int rDimension) :
        Interpolation3D::Interpolation3D(rDofType, rTypeOrder, rDimension)
{
    Initialize();
}

NuTo::eIntegrationType NuTo::Interpolation3DTetrahedron::GetStandardIntegrationType() const
{
    switch (mTypeOrder)
    {
    case NuTo::Interpolation::eTypeOrder::EQUIDISTANT1:
        return NuTo::eIntegrationType::IntegrationType3D4NGauss1Ip;
    case NuTo::Interpolation::eTypeOrder::EQUIDISTANT2:
        return NuTo::eIntegrationType::IntegrationType3D4NGauss4Ip;
    default:
        throw MechanicsException(__PRETTY_FUNCTION__, "Interpolation for exact integration of " + Interpolation::TypeOrderToString(mTypeOrder) + " not implemented");
    }
}

Eigen::VectorXd NuTo::Interpolation3DTetrahedron::CalculateShapeFunctions(const Eigen::VectorXd& rCoordinates) const
{
    switch (mTypeOrder)
    {
    case NuTo::Interpolation::eTypeOrder::EQUIDISTANT1:
        return ShapeFunctions3D::ShapeFunctionsTetrahedronOrder1(rCoordinates);
    case NuTo::Interpolation::eTypeOrder::EQUIDISTANT2:
        return ShapeFunctions3D::ShapeFunctionsTetrahedronOrder2(rCoordinates);
    default:
        throw MechanicsException(__PRETTY_FUNCTION__, "Interpolation order for " + Interpolation::TypeOrderToString(mTypeOrder) + " not implemented");
    }
}

Eigen::MatrixXd NuTo::Interpolation3DTetrahedron::CalculateDerivativeShapeFunctionsNatural(const Eigen::VectorXd& rCoordinates) const
{
    switch (mTypeOrder)
    {
    case NuTo::Interpolation::eTypeOrder::EQUIDISTANT1:
        return ShapeFunctions3D::DerivativeShapeFunctionsTetrahedronOrder1(rCoordinates);
    case NuTo::Interpolation::eTypeOrder::EQUIDISTANT2:
        return ShapeFunctions3D::DerivativeShapeFunctionsTetrahedronOrder2(rCoordinates);
    default:
        throw MechanicsException(__PRETTY_FUNCTION__, "Interpolation order for " + Interpolation::TypeOrderToString(mTypeOrder) + " not implemented");
    }
}

Eigen::VectorXd NuTo::Interpolation3DTetrahedron::CalculateNaturalNodeCoordinates(int rNodeIndexDof) const
{
    switch (mTypeOrder)
    {
    case NuTo::Interpolation::eTypeOrder::EQUIDISTANT1:
        return ShapeFunctions3D::NodeCoordinatesTetrahedronOrder1(rNodeIndexDof);
    case NuTo::Interpolation::eTypeOrder::EQUIDISTANT2:
        return ShapeFunctions3D::NodeCoordinatesTetrahedronOrder2(rNodeIndexDof);
    default:
        throw MechanicsException(__PRETTY_FUNCTION__, "Node arrangement for " + Interpolation::TypeOrderToString(mTypeOrder) + " not implemented");
    }
}

int NuTo::Interpolation3DTetrahedron::CalculateNumNodes() const
{
    switch (mTypeOrder)
    {
    case NuTo::Interpolation::eTypeOrder::EQUIDISTANT1:
        return 4;
    case NuTo::Interpolation::eTypeOrder::EQUIDISTANT2:
        return 10;
    default:
        throw MechanicsException(__PRETTY_FUNCTION__, "Interpolation type and order " + Interpolation::TypeOrderToString(mTypeOrder) + " not implemented");
    }
}

Eigen::VectorXd NuTo::Interpolation3DTetrahedron::CalculateNaturalSurfaceCoordinates(const Eigen::VectorXd& rNaturalSurfaceCoordinates, int rSurface) const
{
    assert(rNaturalSurfaceCoordinates.rows() == 2);
    double alpha = rNaturalSurfaceCoordinates(0);
    double beta = rNaturalSurfaceCoordinates(1);

    switch (rSurface)
    {
    case 0:
        return Eigen::Vector3d(beta, alpha, 0);
    case 1:
        return Eigen::Vector3d(0, beta, alpha);
    case 2:
        return Eigen::Vector3d(alpha, 0, beta);
    case 3:
        return Eigen::Vector3d(1 - alpha - beta, alpha, beta);
    default:
        throw MechanicsException(__PRETTY_FUNCTION__, "TETRAHEDRON3D has exactly four surfaces, 0 to 3. You tried to access " + std::to_string(rSurface) + ".");
    }
}

Eigen::MatrixXd NuTo::Interpolation3DTetrahedron::CalculateDerivativeNaturalSurfaceCoordinates(const Eigen::VectorXd& rNaturalSurfaceCoordinates, int rSurface) const
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
        dXidAlpha(1, 1) = 1.;
        dXidAlpha(2, 0) = 1.;
        break;
    case 2:
        dXidAlpha(0, 0) = 1.;
        dXidAlpha(2, 1) = 1.;
        break;
    case 3:
        dXidAlpha(0, 0) = -1.;
        dXidAlpha(0, 1) = -1.;
        dXidAlpha(1, 0) = 1.;
        dXidAlpha(2, 1) = 1.;
        break;
    default:
        throw MechanicsException(__PRETTY_FUNCTION__, "TETRAHEDRON3D has exactly four surfaces, 0 to 3. You tried to access " + std::to_string(rSurface) + ".");
    }
    return dXidAlpha;
}

std::vector<Eigen::VectorXd> NuTo::Interpolation3DTetrahedron::GetSurfaceEdgesCoordinates(int rSurface) const
{
    int numNodes = 3;
    // returns exactly three nodes, (0,0).T; (1,0).T; and (0,1).T
    Eigen::Matrix<double, 2, 3> alpha = Eigen::Matrix<double, 2, 3>::Zero(); // row1 = alpha, row2 = beta
    alpha(0, 1) = 1;
    alpha(1, 2) = 1;

    std::vector<Eigen::VectorXd> surfaceEdgeCoordinates(numNodes);
    for (int i = 0; i < numNodes; ++i)
    {
        Eigen::VectorXd naturalSurfaceCoordinate = alpha.col(i);
        surfaceEdgeCoordinates[i] = CalculateNaturalSurfaceCoordinates(naturalSurfaceCoordinate, rSurface);
    }
    return surfaceEdgeCoordinates;
}

#ifdef ENABLE_SERIALIZATION
BOOST_CLASS_EXPORT_IMPLEMENT(NuTo::Interpolation3DTetrahedron)
#endif

