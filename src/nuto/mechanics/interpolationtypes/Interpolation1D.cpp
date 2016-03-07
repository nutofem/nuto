/*
 * Interpolation1D.cpp
 *
 *  Created on: 3 Jun 2015
 *      Author: ttitsche
 */

#include "nuto/mechanics/interpolationtypes/Interpolation1D.h"

NuTo::Interpolation1D::Interpolation1D(NuTo::Node::eAttributes rDofType, NuTo::Interpolation::eTypeOrder rTypeOrder, int rDimension) :
        InterpolationBase::InterpolationBase(rDofType, rTypeOrder, rDimension)
{

}

const std::vector<Eigen::VectorXd> NuTo::Interpolation1D::GetSurfaceEdgesCoordinates(int rSurface) const
{
    Eigen::VectorXd dummy; // has no influence in 1D
    return std::vector<Eigen::VectorXd>(1, this->CalculateNaturalSurfaceCoordinates(dummy, rSurface));
}

int NuTo::Interpolation1D::GetNumDofsPerNode() const
{
    switch (mDofType)
    {
    case NuTo::Node::COORDINATES:
        return mDimension;
    case NuTo::Node::DISPLACEMENTS:
        return mDimension;
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

