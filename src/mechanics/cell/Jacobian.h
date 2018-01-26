#pragma once

#include <Eigen/Dense> // for determinant
#include "mechanics/interpolation/TypeDefs.h"

namespace NuTo
{

class Jacobian
{
public:
    using Dynamic3by3 = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::ColMajor, 3, 3>;

    Jacobian(const NodeValues& nodeValues, const DerivativeShapeFunctionsNatural& derivativeShapeFunctions,
             int globalDimension)
    {
        const int interpolationDimension = derivativeShapeFunctions.cols();
        // case 1: global dimension ( node dimension ) matches the interpolation dimension.
        if (interpolationDimension == globalDimension)
        {
            switch (globalDimension)
            {
            case 1:
            {
                Eigen::Matrix<double, 1, 1> jacobian = CalculateFixedSize<1,1>(nodeValues, derivativeShapeFunctions);
                mJacobian = jacobian;
                mInvJacobian = jacobian.inverse();
                mDetJacobian = jacobian.determinant();
                break;
            }
            case 2:
            {
                Eigen::Matrix<double, 2, 2> jacobian = CalculateFixedSize<2,2>(nodeValues, derivativeShapeFunctions);
                mJacobian = jacobian;
                mInvJacobian = jacobian.inverse();
                mDetJacobian = jacobian.determinant();
                break;
            }
            case 3:
            {
                Eigen::Matrix<double, 3, 3> jacobian = CalculateFixedSize<3,3>(nodeValues, derivativeShapeFunctions);
                mJacobian = jacobian;
                mInvJacobian = jacobian.inverse();
                mDetJacobian = jacobian.determinant();
                break;
            }
            }
        }
        else
        {
            assert(globalDimension == interpolationDimension + 1);
            switch (globalDimension)
            {
            case 2:
            {
                const int numRows = derivativeShapeFunctions.rows();
                NodeValues nodeCopy = nodeValues;
                // see CalculateFixedSize(...) for the Eigen::Map functionality
                auto nodeBlockCoordinates =
                        Eigen::Map<Eigen::Matrix<double, 2, Eigen::Dynamic>>(nodeCopy.data(), 2, numRows);
                mDetJacobian = (nodeBlockCoordinates * derivativeShapeFunctions).norm();
                break;
            }
            case 3:
            {
                const int numRows = derivativeShapeFunctions.rows();
                NodeValues nodeCopy = nodeValues;
                // see CalculateFixedSize(...) for the Eigen::Map functionality
                auto nodeBlockCoordinatesMap =
                        Eigen::Map<Eigen::Matrix<double, 3, Eigen::Dynamic>>(nodeCopy.data(), 3, numRows);

                Eigen::MatrixXd nodeBlockCoordinates = nodeBlockCoordinatesMap;

                Eigen::Matrix<double, 3, 2> jacobian = nodeBlockCoordinates * derivativeShapeFunctions;
                Eigen::Vector3d derivativeXi = jacobian.col(0);
                Eigen::Vector3d derivativeMu = jacobian.col(1);

                mDetJacobian = (derivativeXi.cross(derivativeMu)).norm();
                break;
            }
            default:
            {
                throw;
            }
            }
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

    const Dynamic3by3& Get() const
    {
        return mJacobian;
    }

    double Det() const
    {
        return mDetJacobian;
    }

private:
    template <int TSpaceDim, int TLocalDim>
    Eigen::Matrix<double, TSpaceDim, TLocalDim>
    CalculateFixedSize(NuTo::NodeValues nodeValues,
                       const NuTo::DerivativeShapeFunctionsNatural& derivativeShapeFunctions)
    {
        const int numShapes = derivativeShapeFunctions.rows();
        auto nodeBlockCoordinates =
                Eigen::Map<Eigen::Matrix<double, TSpaceDim, Eigen::Dynamic>>(nodeValues.data(), TSpaceDim, numShapes);
        // This mapping converts a vector
        // (x0 y0 z0 x1 y1 z1 x2 y2 z3 ...)^T
        // to the matrix
        // x0  x1  x2  x3 ...
        // y0  y1  y2  y3 ...
        // z0  z1  z2  z3 ...
        Eigen::Matrix<double, TSpaceDim, TLocalDim> jacobian = nodeBlockCoordinates * derivativeShapeFunctions;
        return jacobian;
    }

    Dynamic3by3 mJacobian;
    Dynamic3by3 mInvJacobian;
    double mDetJacobian;
};
} /* NuTo */
