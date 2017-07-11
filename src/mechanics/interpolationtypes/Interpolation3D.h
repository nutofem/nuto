/*
 * Interpolation3D.h
 *
 *  Created on: 19 May 2015
 *      Author: ttitsche
 */

#pragma once

#include "mechanics/interpolationtypes/InterpolationBaseFEM.h"

namespace NuTo
{

class Interpolation3D : public InterpolationBaseFEM
{
public:
    Interpolation3D(NuTo::Node::eDof rDofType, NuTo::Interpolation::eTypeOrder rTypeOrder, int rDimension);

    virtual int GetLocalDimension() const override
    {
        return 3;
    }
};

} /* namespace NuTo */
