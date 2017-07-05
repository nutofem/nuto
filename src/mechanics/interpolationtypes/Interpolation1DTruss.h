/*
 * Interpolation1DTruss.h
 *
 *  Created on: 8 May 2015
 *      Author: ttitsche
 */

#pragma once

#include "mechanics/interpolationtypes/Interpolation1D.h"

namespace NuTo
{

/**
@brief 2D quadrilateral element with the following natural coordinate system and its surface parametrization
@image html Truss1D.png
**/
class Interpolation1DTruss : public Interpolation1D
{

public:
    Interpolation1DTruss(NuTo::Node::eDof rDofType, NuTo::Interpolation::eTypeOrder rTypeOrder, int rDimension);

    eIntegrationType GetStandardIntegrationType() const override;

    Eigen::VectorXd CalculateNaturalSurfaceCoordinates(const Eigen::VectorXd& rNaturalSurfaceCoordinates,
                                                       int rSurface) const override;

    Eigen::MatrixXd CalculateDerivativeNaturalSurfaceCoordinates(const Eigen::VectorXd& rNaturalSurfaceCoordinates,
                                                                 int rSurface) const override;


    inline int GetNumSurfaces() const override
    {
        return 2;
    }

protected:
    Eigen::VectorXd CalculateNaturalNodeCoordinates(int rNodeIndexDof) const override;

    Eigen::VectorXd CalculateShapeFunctions(const Eigen::VectorXd& rCoordinates) const override;

    Eigen::MatrixXd CalculateDerivativeShapeFunctionsNatural(const Eigen::VectorXd& rCoordinates) const override;

private:
    int CalculateNumNodes() const override;
};

} /* namespace NuTo */
