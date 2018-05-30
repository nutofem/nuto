#pragma once

#include <vector>
#include "nuto/mechanics/nodes/DofNode.h"
#include "nuto/mechanics/elements/DofElementInterface.h"
#include "nuto/mechanics/iga/Nurbs.h"
#include "nuto/mechanics/cell/Matrix.h"

namespace NuTo
{
template <int TDimParameter>
class ElementIga : public DofElementInterface
{
public:
    ElementIga(const std::array<int, TDimParameter>& knotIDs, const Nurbs<TDimParameter>& NurbsGeometry)
        : mKnotIDs(knotIDs)
        , mNurbsGeometry(NurbsGeometry)
    {
    }

    //! @brief transforms unit interval [-1, 1] to the interval [firstKnotCoordinate, secondKnotCoordinate]
    Eigen::Matrix<double, TDimParameter, 1> Transformation(Eigen::VectorXd ipCoords) const
    {
        Eigen::Matrix<double, TDimParameter, 1> coordinateTransformed;
        std::array<Eigen::Vector2d, TDimParameter> knots = NurbsGeometry().GetKnotVectorElement(mKnotIDs);
        for (int i = 0; i < TDimParameter; i++)
        {
            coordinateTransformed[i] = (knots[i](0) + 0.5 * (ipCoords(i) + 1) * (knots[i](1) - knots[i](0)));
        }
        return coordinateTransformed;
    }

    //! @brief extracts all node values of this element
    Eigen::VectorXd ExtractNodeValues(int instance = 0) const override
    {
        return NurbsGeometry().GetControlPointsElement(mKnotIDs, instance);
    }

    Eigen::MatrixXd GetNMatrix(NaturalCoords ipCoords) const override
    {
        return NuTo::Matrix::N(GetShapeFunctions(ipCoords), GetNumNodes(), GetDofDimension());
    }

    Eigen::VectorXd GetShapeFunctions(NaturalCoords ipCoords) const override
    {
        return NurbsGeometry().BasisFunctionsAndDerivativesRational(0, Transformation(ipCoords));
    }

    Eigen::MatrixXd GetDerivativeShapeFunctions(NaturalCoords ipCoords) const override
    {
        return NurbsGeometry().BasisFunctionsAndDerivativesRational(1, Transformation(ipCoords));
    }

    int GetDofDimension() const override
    {
        return NurbsGeometry().GetDimension();
    }

    Eigen::VectorXi GetDofNumbering() const override
    {
        throw NuTo::Exception(__PRETTY_FUNCTION__, "I honestly have no idea. Sorry.");
    }

    int GetNumNodes() const override
    {
        return NurbsGeometry().GetNumControlPointsElement();
    }

private:
    const Nurbs<TDimParameter>& NurbsGeometry() const
    {
        return mNurbsGeometry;
    }

    std::array<int, TDimParameter> mKnotIDs;
    std::reference_wrapper<const Nurbs<TDimParameter>> mNurbsGeometry;
};
} /* NuTo */
