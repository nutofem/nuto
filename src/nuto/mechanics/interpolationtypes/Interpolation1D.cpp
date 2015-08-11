/*
 * Interpolation1D.cpp
 *
 *  Created on: 3 Jun 2015
 *      Author: ttitsche
 */

#include "nuto/mechanics/interpolationtypes/Interpolation1D.h"

NuTo::Interpolation1D::Interpolation1D(const StructureBase* rStructure, NuTo::Node::eAttributes rDofType, NuTo::Interpolation::eTypeOrder rTypeOrder) :
        InterpolationBase::InterpolationBase(rStructure, rDofType, rTypeOrder)
{

}

const std::vector<Eigen::VectorXd> NuTo::Interpolation1D::GetSurfaceEdgesCoordinates(int rSurface) const
{
    Eigen::VectorXd dummy; // has no influence in 1D
    return std::vector<Eigen::VectorXd>(1, CalculateNaturalSurfaceCoordinates(dummy, rSurface));
}

int NuTo::Interpolation1D::GetNumDofsPerNode() const
{
    switch (mDofType)
    {
    case NuTo::Node::COORDINATES:
        return mStructure->GetDimension();
    case NuTo::Node::DISPLACEMENTS:
        return mStructure->GetDimension();
    case NuTo::Node::TEMPERATURES:
        return 1;
    case NuTo::Node::NONLOCALEQSTRAIN:
        return 1;
    case NuTo::Node::NONLOCALEQPLASTICSTRAIN:
        return 2;
    case NuTo::Node::RELATIVEHUMIDITY:
        return 1;
    case NuTo::Node::WATERVOLUMEFRACTION:
        return 1;
    default:
        throw NuTo::MechanicsException("[NuTo::Interpolation1D::GetNumDofsPerNode] dof type not found.");
    }
}

