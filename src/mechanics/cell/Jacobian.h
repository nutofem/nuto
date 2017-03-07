#pragma once

#include <eigen3/Eigen/Dense>
#include "mechanics/interpolation/TypeDefs.h"

namespace NuTo
{

template <int TDim>
class Jacobian
{
public:
    Jacobian(const NodeValues& rNodeValues, const DerivativeShapeFunctionsNatural& rDerivativeShapeFunctions)
    {
        const int numRows = rDerivativeShapeFunctions.rows();
        Eigen::MatrixXd nodeBlockCoordinates(TDim, numRows);
        // convert the coordinates to a block structure
        // x0  x1  x2  x3 ...
        // y0  y1  y2  y3 ...
        // z0  z1  z2  z3 ...
        for (int i                      = 0; i < numRows; ++i)
            nodeBlockCoordinates.col(i) = rNodeValues.block<TDim, 1>(TDim * i, 0);

        mJacobian    = nodeBlockCoordinates.lazyProduct(rDerivativeShapeFunctions);
        mInvJacobian = mJacobian.inverse();
        mDetJacobian = mJacobian.determinant();
    }

    DerivativeShapeFunctionsGlobal
    TransformDerivativeShapeFunctions(const DerivativeShapeFunctionsNatural& rGlobal) const
    {
        return rGlobal * Inv();
    }

    const Eigen::Matrix<double, TDim, TDim>& Inv() const
    {
        return mInvJacobian;
    }

    double Det() const
    {
        return mDetJacobian;
    }

private:
    Eigen::Matrix<double, TDim, TDim> mJacobian;
    Eigen::Matrix<double, TDim, TDim> mInvJacobian;
    double mDetJacobian;
};
} /* NuTo */
