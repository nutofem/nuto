#pragma once

#include <Eigen/Dense> // for determinant
#include "nuto/mechanics/interpolation/TypeDefs.h"
#include "nuto/base/Exception.h"

namespace NuTo
{

class Jacobian
{
public:
    using Dynamic3by3 = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::ColMajor, 3, 3>;
    using Dynamic3by1 = Eigen::Matrix<double, Eigen::Dynamic, 1, Eigen::ColMajor, 3, 1>;

    Jacobian(const Eigen::VectorXd& nodeValues, const Eigen::MatrixXd& derivativeShapeFunctions, int globalDimension)
    {
        const int interpolationDimension = derivativeShapeFunctions.cols();
        // case 1: global dimension ( node dimension ) matches the interpolation dimension.
        if (interpolationDimension == globalDimension)
        {
            switch (globalDimension)
            {
            case 1:
            {
                Eigen::Matrix<double, 1, 1> jacobian = CalculateFixedSize<1, 1>(nodeValues, derivativeShapeFunctions);
                mJacobian = jacobian;
                mInvJacobian = jacobian.inverse();
                mDetJacobian = jacobian.determinant();
                break;
            }
            case 2:
            {
                Eigen::Matrix<double, 2, 2> jacobian = CalculateFixedSize<2, 2>(nodeValues, derivativeShapeFunctions);
                mJacobian = jacobian;
                mInvJacobian = jacobian.inverse();
                mDetJacobian = jacobian.determinant();
                break;
            }
            case 3:
            {
                Eigen::Matrix<double, 3, 3> jacobian = CalculateFixedSize<3, 3>(nodeValues, derivativeShapeFunctions);
                mJacobian = jacobian;
                mInvJacobian = jacobian.inverse();
                mDetJacobian = jacobian.determinant();
                break;
            }
            }
        }
        else if (globalDimension == interpolationDimension + 1)
        {
            switch (globalDimension)
            {
            case 2:
            {
                Eigen::Matrix<double, 2, 1> jacobian = CalculateFixedSize<2, 1>(nodeValues, derivativeShapeFunctions);
                Eigen::Matrix<double, 2, 2> extendedJacobian;
                mNormal = Eigen::Vector2d(jacobian(1, 0), -jacobian(0, 0)).normalized();
                extendedJacobian.col(0) = jacobian;
                extendedJacobian.col(1) = mNormal;
                mJacobian = jacobian;
                mInvJacobian = extendedJacobian.inverse().block<1, 2>(0, 0);
                mDetJacobian = -extendedJacobian.determinant();
                break;
            }
            case 3:
            {
                Eigen::Matrix<double, 3, 2> jacobian = CalculateFixedSize<3, 2>(nodeValues, derivativeShapeFunctions);
                Eigen::Matrix<double, 3, 3> extendedJacobian;
                mNormal = (jacobian.col(0)).cross(jacobian.col(1)).normalized();
                extendedJacobian.col(0) = jacobian.col(0);
                extendedJacobian.col(1) = jacobian.col(1);
                extendedJacobian.col(2) = mNormal;
                mJacobian = jacobian;
                mInvJacobian = extendedJacobian.inverse().block<2, 3>(0, 0);
                ;
                mDetJacobian = extendedJacobian.determinant();
                break;
            }
            default:
            {
                throw Exception(__PRETTY_FUNCTION__,
                                "Global dimension should be 2 or 3 for codimension 1 objects, got: " +
                                        std::to_string(globalDimension) + ".");
            }
            }
        }
        else
        {
            throw Exception(__PRETTY_FUNCTION__,
                            "Jacobian only implemented for interpolationDimension equal to globaldimension or one "
                            "less, got: \n"
                            "InterpolationDimension: " +
                                    std::to_string(interpolationDimension) + "\n" +
                                    "GlobalDimension: " + std::to_string(globalDimension));
        }
    }

    Eigen::MatrixXd TransformDerivativeShapeFunctions(const Eigen::MatrixXd& global) const
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

    const Dynamic3by1& Normal() const
    {
        if (mJacobian.cols() != (mJacobian.rows() - 1))
        {
            throw Exception(__PRETTY_FUNCTION__,
                            "Normal not available for elements with SpaceDimension - LocalDimension != 1");
        }
        return mNormal;
    }

    double Det() const
    {
        return mDetJacobian;
    }

private:
    template <int TSpaceDim, int TLocalDim>
    Eigen::Matrix<double, TSpaceDim, TLocalDim> CalculateFixedSize(Eigen::VectorXd nodeValues,
                                                                   const Eigen::MatrixXd& derivativeShapeFunctions)
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

    Dynamic3by1 mNormal;
    Dynamic3by3 mJacobian;
    Dynamic3by3 mInvJacobian;
    double mDetJacobian;
};
} /* NuTo */
