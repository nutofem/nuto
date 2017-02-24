#pragma once
#include <eigen3/Eigen/Core>
#include "mechanics/interpolation/Interpolation.h"
#include "mechanics/elements/ElementShapeFunctions.h"

namespace NuTo
{
class InterpolationTriangle : public Interpolation
{
public:
    InterpolationTriangle(eInterpolation rType, int rOrder)
        : Interpolation(rType, rOrder)
    {
    }

    Eigen::VectorXd GetShapeFunctions(const Eigen::VectorXd& rIPCoords) const override
    {
        if (mType == eInterpolation::GAUSS && mOrder == 1)
            return ShapeFunctions2D::ShapeFunctionsTriangleOrder1(rIPCoords);
        throw;
    }

    Eigen::MatrixXd GetDerivativeShapeFunctions(const Eigen::VectorXd& rIPCoords) const override
    {
        if (mType == eInterpolation::GAUSS && mOrder == 1)
            return ShapeFunctions2D::DerivativeShapeFunctionsTriangleOrder1(rIPCoords);
        throw;
    }

    Eigen::VectorXd GetLocalCoords(int rNodeId) const override
    {
        if (mType == eInterpolation::GAUSS && mOrder == 1)
            return ShapeFunctions2D::NodeCoordinatesTriangleOrder1(rNodeId);
        throw;
    }
    int GetNumNodes() const override
    {
        if (mType == eInterpolation::GAUSS && mOrder == 1)
            return 3;
        throw;
    }


private:
};
} /* NuTo */
