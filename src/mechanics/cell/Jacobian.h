#pragma once

#include <eigen3/Eigen/Dense> // for determinant
#include "mechanics/interpolation/TypeDefs.h"

namespace NuTo
{

class Jacobian
{
public:
    using Dynamic3by3 = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::ColMajor, 3, 3>;

    Jacobian(const NodeValues& nodeValues, const DerivativeShapeFunctionsNatural& derivativeShapeFunctions)
    {
        const int dim = nodeValues.rows() / derivativeShapeFunctions.rows();
        switch (dim)
        {
        case 1:
            std::tie(mInvJacobian, mDetJacobian) = CalculateFixedSize<1>(nodeValues, derivativeShapeFunctions);
            break;
        case 2:
            std::tie(mInvJacobian, mDetJacobian) = CalculateFixedSize<2>(nodeValues, derivativeShapeFunctions);
            break;
        case 3:
            std::tie(mInvJacobian, mDetJacobian) = CalculateFixedSize<3>(nodeValues, derivativeShapeFunctions);
            break;
        }
    }

    DerivativeShapeFunctionsGlobal
    TransformDerivativeShapeFunctions(const DerivativeShapeFunctionsNatural& global) const
    {
        return global * Inv();
    }

    const Dynamic3by3& Inv() const
    {
        return mInvJacobian;
    }

    double Det() const
    {
        return mDetJacobian;
    }

private:
    template <int TDim>
    std::pair<Eigen::Matrix<double, TDim, TDim>, double>
    CalculateFixedSize(NuTo::NodeValues nodeValues,
                       const NuTo::DerivativeShapeFunctionsNatural& derivativeShapeFunctions)
    {
        const int numRows = derivativeShapeFunctions.rows();
        auto nodeBlockCoordinates =
                Eigen::Map<Eigen::Matrix<double, TDim, Eigen::Dynamic>>(nodeValues.data(), TDim, numRows);
        // This mapping converts a vector
        // (x0 y0 z0 x1 y1 z1 x2 y2 z3 ...)^T
        // to the matrix
        // x0  x1  x2  x3 ...
        // y0  y1  y2  y3 ...
        // z0  z1  z2  z3 ...
        Eigen::Matrix<double, TDim, TDim> jacobian = nodeBlockCoordinates * derivativeShapeFunctions;
        return std::make_pair(jacobian.inverse(), jacobian.determinant());
    }

    Dynamic3by3 mInvJacobian;
    double mDetJacobian;
};
} /* NuTo */
