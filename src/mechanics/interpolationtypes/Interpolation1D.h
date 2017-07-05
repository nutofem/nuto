/*
 * Interpolation1D.h
 *
 *  Created on: 3 Jun 2015
 *      Author: ttitsche
 */

#pragma once

#include "mechanics/interpolationtypes/InterpolationBaseFEM.h"

namespace NuTo
{

class Interpolation1D: public InterpolationBaseFEM
{

public:

    Interpolation1D(NuTo::Node::eDof rDofType, NuTo::Interpolation::eTypeOrder rTypeOrder, int rDimension);

    std::vector<Eigen::VectorXd> GetSurfaceEdgesCoordinates(int rSurface) const override;

    virtual int GetLocalDimension() const override
    {
        return 1;
    }

};

} /* namespace NuTo */


