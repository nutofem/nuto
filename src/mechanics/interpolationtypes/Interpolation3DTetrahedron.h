/*
 * Interpolation3DTetrahedron.h
 *
 *  Created on: 19 May 2015
 *      Author: ttitsche
 */

#pragma once

#include "mechanics/interpolationtypes/Interpolation3D.h"

namespace NuTo
{
/**
@brief 3D tetrahedral element with the following natural coordinate system and its surface parametrization
@image html Tetrahedron3D.png
**/
class Interpolation3DTetrahedron: public Interpolation3D
{

public:
    Interpolation3DTetrahedron(NuTo::Node::eDof rDofType, NuTo::Interpolation::eTypeOrder rTypeOrder, int rDimension);

    eIntegrationType GetStandardIntegrationType() const override;

    Eigen::VectorXd CalculateNaturalSurfaceCoordinates(const Eigen::VectorXd& rNaturalSurfaceCoordinates, int rSurface) const override;

    Eigen::MatrixXd CalculateDerivativeNaturalSurfaceCoordinates(const Eigen::VectorXd& rNaturalSurfaceCoordinates, int rSurface) const override;

    inline int GetNumSurfaces() const override
    {
        return 4;
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

