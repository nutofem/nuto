#pragma once
#include <eigen3/Eigen/Core>
#include "mechanics/interpolation/InterpolationEnum.h"

namespace NuTo
{
class Interpolation
{
public:
    Interpolation(eInterpolation rType, int rOrder)
        : mType(rType)
        , mOrder(rOrder)
    {
    }
    virtual ~Interpolation() = default;

    Interpolation(const Interpolation&) = default;
    Interpolation(Interpolation&&)      = default;
    Interpolation& operator=(const Interpolation&) = default;
    Interpolation& operator=(Interpolation&&) = default;

    virtual Eigen::VectorXd GetShapeFunctions(const Eigen::VectorXd& rIPCoords) const = 0;
    virtual Eigen::MatrixXd GetDerivativeShapeFunctions(const Eigen::VectorXd& rIPCoords) const = 0;
    virtual Eigen::VectorXd GetLocalCoords(int rNodeId) const            = 0;
    virtual int GetNumNodes() const                                      = 0;

protected:
    eInterpolation mType;
    int mOrder;
};
} /* NuTo */
