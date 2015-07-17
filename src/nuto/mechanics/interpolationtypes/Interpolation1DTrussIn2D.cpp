/*
 * Interpolation1DTrussIn2D.cpp
 *
 *  Created on: 10 July 2015
 *      Author: phuschke
 */

#include "nuto/mechanics/interpolationtypes/Interpolation1DTrussIn2D.h"

NuTo::Interpolation1DTrussIn2D::Interpolation1DTrussIn2D(NuTo::Node::eAttributes rDofType, NuTo::Interpolation::eTypeOrder rTypeOrder) :
        NuTo::Interpolation1DTruss::Interpolation1DTruss(rDofType, rTypeOrder)
{
    Initialize();
}

int NuTo::Interpolation1DTrussIn2D::GetNumDofsPerNode() const
{
    return 2;
}
