/*
 * Interpolation3DPrism.h
 *
 *  Created on: 10 January 2017
 *      Author: Thomas Titscher
 */

#pragma once

#include "mechanics/interpolationtypes/Interpolation3D.h"

namespace NuTo
{
/**
@brief 3D prism element with the following natural coordinate system and its surface parametrization
@image html Prism3D.png
**/
class Interpolation3DPrism : public Interpolation3D
{
public:
    Interpolation3DPrism(NuTo::Node::eDof rDofType, NuTo::Interpolation::eTypeOrder rTypeOrder, int rDimension);

    eIntegrationType GetStandardIntegrationType() const override;

    Eigen::VectorXd CalculateNaturalSurfaceCoordinates(const Eigen::VectorXd& rNaturalSurfaceCoordinates,
                                                       int rSurface) const override;

    Eigen::MatrixXd CalculateDerivativeNaturalSurfaceCoordinates(const Eigen::VectorXd& rNaturalSurfaceCoordinates,
                                                                 int rSurface) const override;

    int GetNumSurfaces() const override
    {
        return 5;
    }

    std::vector<Eigen::VectorXd> GetSurfaceEdgesCoordinates(int rSurface) const override;

protected:
    Eigen::VectorXd CalculateNaturalNodeCoordinates(int rNodeIndexDof) const override;

    Eigen::VectorXd CalculateShapeFunctions(const Eigen::VectorXd& rCoordinates) const override;

    Eigen::MatrixXd CalculateDerivativeShapeFunctionsNatural(const Eigen::VectorXd& rCoordinates) const override;

private:
    int CalculateNumNodes() const override;
};

} /* namespace NuTo */
