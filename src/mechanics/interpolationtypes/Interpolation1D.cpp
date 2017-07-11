/*
 * Interpolation1D.cpp
 *
 *  Created on: 3 Jun 2015
 *      Author: ttitsche
 */

#include "mechanics/MechanicsException.h"
#include "mechanics/interpolationtypes/Interpolation1D.h"

NuTo::Interpolation1D::Interpolation1D(NuTo::Node::eDof rDofType, NuTo::Interpolation::eTypeOrder rTypeOrder,
                                       int rDimension)
    : InterpolationBaseFEM::InterpolationBaseFEM(rDofType, rTypeOrder, rDimension)
{
}

std::vector<Eigen::VectorXd> NuTo::Interpolation1D::GetSurfaceEdgesCoordinates(int rSurface) const
{
    Eigen::VectorXd dummy; // has no influence in 1D
    return std::vector<Eigen::VectorXd>(1, this->CalculateNaturalSurfaceCoordinates(dummy, rSurface));
}
