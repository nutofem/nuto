/*
 * InterpolationType3D.cpp
 *
 *  Created on: 19 May 2015
 *      Author: ttitsche
 */

#include "base/Exception.h"
#include "mechanics/interpolationtypes/Interpolation3D.h"

NuTo::Interpolation3D::Interpolation3D(NuTo::Node::eDof rDofType, NuTo::Interpolation::eTypeOrder rTypeOrder,
                                       int rDimension)
    : InterpolationBaseFEM::InterpolationBaseFEM(rDofType, rTypeOrder, rDimension)
{
}
