/*
 * InterpolationType2D.cpp
 *
 *  Created on: 23 Mar 2015
 *      Author: ttitsche
 */

#include "nuto/mechanics/interpolationtypes/Interpolation2D.h"

NuTo::Interpolation2D::Interpolation2D(NuTo::Node::eAttributes rDofType, NuTo::Interpolation::eTypeOrder rTypeOrder) :
        InterpolationBase::InterpolationBase(rDofType, rTypeOrder)
{

}

const std::vector<Eigen::VectorXd> NuTo::Interpolation2D::GetSurfaceEdgesCoordinates(int rSurface) const
{
    int numNodes = 2;
    // returns exactly two nodes, one at alpha = -1 and one at alpha = 1 of the required surface
    Eigen::Vector2d alpha(-1, 1);
    std::vector<Eigen::VectorXd> surfaceEdgeCoordinates(numNodes);
    for (int i = 0; i < numNodes; ++i)
    {
        Eigen::VectorXd naturalSurfaceCoordinate = alpha.row(i);
        surfaceEdgeCoordinates[i] = CalculateNaturalSurfaceCoordinates(naturalSurfaceCoordinate, rSurface);
    }
    return surfaceEdgeCoordinates;
}

int NuTo::Interpolation2D::GetNumDofsPerNode() const
{
    switch (mDofType)
    {
    case NuTo::Node::COORDINATES:
        return 2;
    case NuTo::Node::DISPLACEMENTS:
        return 2;
    case NuTo::Node::TEMPERATURES:
        return 1;
    case NuTo::Node::NONLOCALEQSTRAIN:
        return 1;
    case NuTo::Node::NONLOCALEQPLASTICSTRAIN:
        return 2;
    default:
        throw NuTo::MechanicsException("[NuTo::Interpolation2D::GetNumDofsPerNode] dof type not found.");
    }
}
