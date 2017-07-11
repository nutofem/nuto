/*
 * Interpolation2D.h
 *
 *  Created on: 23 Mar 2015
 *      Author: ttitsche
 */

#pragma once

#include "mechanics/interpolationtypes/InterpolationBaseFEM.h"

namespace NuTo
{

class Interpolation2D : public InterpolationBaseFEM
{

public:
    Interpolation2D(NuTo::Node::eDof rDofType, NuTo::Interpolation::eTypeOrder rTypeOrder, int rDimension);

    std::vector<Eigen::VectorXd> GetSurfaceEdgesCoordinates(int rSurface) const override;

    Eigen::VectorXi GetSurfaceNodeIndices(int rSurface) const override
    {
        throw NuTo::Exception(__PRETTY_FUNCTION__,
                                       "Implemented in NuTo::InterpolationType::GetSurfaceNodeIndices.");
    }

    virtual int GetLocalDimension() const override
    {
        return 2;
    }
};

} /* namespace NuTo */
