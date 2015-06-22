/*
 * Interpolation1D.cpp
 *
 *  Created on: 3 Jun 2015
 *      Author: ttitsche
 */

#include "nuto/mechanics/interpolationtypes/Interpolation1D.h"

NuTo::Interpolation1D::Interpolation1D(NuTo::Node::eAttributes rDofType, NuTo::Interpolation::eTypeOrder rTypeOrder) :
        InterpolationBase::InterpolationBase(rDofType, rTypeOrder)
{

}

const std::vector<Eigen::VectorXd> NuTo::Interpolation1D::GetSurfaceEdgesCoordinates(int rSurface) const
{
    Eigen::VectorXd dummy; // has no influence in 1D
    return std::vector<Eigen::VectorXd>(1, CalculateNaturalSurfaceCoordinates(dummy, rSurface));
}

int NuTo::Interpolation1D::GetNumDofsPerNode() const
{
    if (mDofType == Node::NONLOCALEQPLASTICSTRAIN)
        return 2;
    return 1;
}
