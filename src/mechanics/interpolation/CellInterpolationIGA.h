#pragma once

#include <vector>
#include "mechanics/nodes/NodeSimple.h"
#include "CellInterpolationBase.h"
#include "mechanics/IGA/NURBS.h"
#include "mechanics/cell/Matrix.h"

namespace NuTo
{
template <int TDimParameter>
class CellInterpolationIGA : public CellInterpolationBase
{
public:
    CellInterpolationIGA(const std::array<int, TDimParameter>& knotIDs, const NURBS<TDimParameter>& NURBSGeometry)
        : mKnotIDs(knotIDs)
        , mNURBSGeometry(NURBSGeometry)
    {
    }

    //! @brief transforms unit interval [-1, 1] to the interval [firstKnotCoordinate, secondKnotCoordinate]
    Eigen::Matrix<double, TDimParameter, 1> Transformation(Eigen::VectorXd ipCoords) const
    {
        Eigen::Matrix<double, TDimParameter, 1> coordinateTransformed;
        std::array<Eigen::Vector2d, TDimParameter> knots = mNURBSGeometry.GetKnotVectorElement(mKnotIDs);
        for (int i = 0; i < TDimParameter; i++)
        {
            coordinateTransformed[i] = (knots[i](0) + 0.5 * (ipCoords(i) + 1) * (knots[i](1) - knots[i](0)));
        }
        return coordinateTransformed;
    }

    //! @brief extracts all node values of this element
    //! @remark virtual to make it testable
    virtual NodeValues ExtractNodeValues() const override
    {
        return mNURBSGeometry.GetControlPointsElement(mKnotIDs);
    }

    NMatrix GetNMatrix(NaturalCoords ipCoords) const override
    {
        return NuTo::Matrix::N(GetShapeFunctions(ipCoords), GetNumNodes(), GetDofDimension());
    }

    Eigen::VectorXd GetShapeFunctions(Eigen::VectorXd ipCoords) const override
    {
        return mNURBSGeometry.BasisFunctionsAndDerivativesRational(0, Transformation(ipCoords));
    }

    Eigen::MatrixXd GetDerivativeShapeFunctions(Eigen::VectorXd ipCoords) const override
    {
        return mNURBSGeometry.BasisFunctionsAndDerivativesRational(1, Transformation(ipCoords));
    }

    int GetDofDimension() const override
    {
        return mNURBSGeometry.GetDimension();
    }

    int GetNumNodes() const override
    {
        return mNURBSGeometry.GetNumControlPointsElement();
    }

private:
    std::array<int, TDimParameter> mKnotIDs;
    const NURBS<TDimParameter>& mNURBSGeometry;
};
} /* NuTo */
