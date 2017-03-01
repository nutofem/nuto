#pragma once
#include <eigen3/Eigen/Core>
#include "mechanics/interpolation/InterpolationEnum.h"

namespace NuTo
{
class Interpolation
{
public:
    Interpolation(eInterpolation rType, int rOrder, int rDofDimension)
        : mType(rType)
        , mOrder(rOrder)
        , mDofDimension(rDofDimension)
    {
    }
    virtual ~Interpolation() = default;

    Interpolation(const Interpolation&) = default;
    Interpolation(Interpolation&&)      = default;
    Interpolation& operator=(const Interpolation&) = default;
    Interpolation& operator=(Interpolation&&) = default;

    int GetDofDimension() const
    {
        return mDofDimension;
    }

    Eigen::MatrixXd GetN(const Eigen::VectorXd& rIPCoords) const;
    Eigen::MatrixXd GetBGradient(const Eigen::VectorXd& rIPCoords) const;
    Eigen::MatrixXd GetBStrain(const Eigen::VectorXd& rIPCoords) const;

    virtual Eigen::VectorXd GetShapeFunctions(const Eigen::VectorXd& rIPCoords) const           = 0;
    virtual Eigen::MatrixXd GetDerivativeShapeFunctions(const Eigen::VectorXd& rIPCoords) const = 0;
    virtual Eigen::VectorXd GetLocalCoords(int rNodeId) const                                   = 0;
    virtual int GetNumNodes() const                                                             = 0;


protected:
    eInterpolation mType;
    int mOrder;
    int mDofDimension;
};
} /* NuTo */
