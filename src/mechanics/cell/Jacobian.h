#pragma once

#include <iostream>
#include <eigen3/Eigen/Dense> // for determinant
#include "mechanics/interpolation/TypeDefs.h"

namespace NuTo
{

class Jacobian
{
public:
    Jacobian(const NodeValues& nodeValues, const DerivativeShapeFunctionsNatural& derivativeShapeFunctions)
    {
        const int numNodes = derivativeShapeFunctions.rows();
        const int dim = derivativeShapeFunctions.cols();

        Eigen::MatrixXd nodeBlockCoordinates(dim, numNodes);
        // convert the coordinates to a block structure
        // x0  x1  x2  x3 ...
        // y0  y1  y2  y3 ...
        // z0  z1  z2  z3 ...
        for (int i = 0; i < numNodes; ++i)
            nodeBlockCoordinates.col(i) = nodeValues.block(dim * i, 0, dim,1);

        mJacobian = nodeBlockCoordinates.lazyProduct(derivativeShapeFunctions);
        mInvJacobian = mJacobian.inverse().eval();
        mDetJacobian = mJacobian.determinant();
    }

    DerivativeShapeFunctionsGlobal
    TransformDerivativeShapeFunctions(const DerivativeShapeFunctionsNatural& global) const
    {
        return global * mInvJacobian;
    }
    const Eigen::MatrixXd& Jac() const
    {
        return mJacobian;
    }

    const Eigen::MatrixXd& Inv() const
    {
        return mInvJacobian;
    }

    double Det() const
    {
        return mDetJacobian;
    }

private:
    Eigen::MatrixXd mJacobian;
    Eigen::MatrixXd mInvJacobian;
    double mDetJacobian;
};
} /* NuTo */
