#pragma once
#include <eigen3/Eigen/Core>
#include "mechanics/interpolation/InterpolationEnum.h"

namespace NuTo
{

typedef Eigen::MatrixXd NMatrix;
typedef Eigen::MatrixXd BMatrixGradient;
typedef Eigen::MatrixXd BMatrixStrain;

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

    NMatrix GetN(const Eigen::VectorXd& rIPCoords) const;
    BMatrixGradient GetBGradient(const Eigen::VectorXd& rIPCoords) const;
    BMatrixStrain GetBStrain(const Eigen::VectorXd& rIPCoords) const;

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
