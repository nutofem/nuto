/*
 * Interpolation2DTriangle.h
 *
 *  Created on: 20 Mar 2015
 *      Author: ttitsche
 */

#pragma once

#include "mechanics/interpolationtypes/Interpolation2D.h"

namespace NuTo
{
/**
@brief 2D triangular element with the following natural coordinate system and its surface parametrization
@image html Triangle2D.png
**/
class Interpolation2DTriangle : public Interpolation2D
{
public:
    Interpolation2DTriangle(NuTo::Node::eDof rDofType, NuTo::Interpolation::eTypeOrder rTypeOrder, int rDimension);

    eIntegrationType GetStandardIntegrationType() const override;

    Eigen::VectorXd CalculateNaturalSurfaceCoordinates(const Eigen::VectorXd& rNaturalSurfaceCoordinates,
                                                       int rSurface) const override;

    Eigen::MatrixXd CalculateDerivativeNaturalSurfaceCoordinates(const Eigen::VectorXd& rNaturalSurfaceCoordinates,
                                                                 int rSurface) const override;

    int GetNumSurfaces() const override
    {
        return 3;
    }

protected:
    Eigen::VectorXd CalculateNaturalNodeCoordinates(int rNodeIndexDof) const override;

    Eigen::VectorXd CalculateShapeFunctions(const Eigen::VectorXd& rCoordinates) const override;

    Eigen::MatrixXd CalculateDerivativeShapeFunctionsNatural(const Eigen::VectorXd& rCoordinates) const override;

private:
    int CalculateNumNodes() const override;
};

} /* namespace NuTo */
