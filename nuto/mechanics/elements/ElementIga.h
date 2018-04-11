#pragma once

#include <vector>
#include "nuto/mechanics/nodes/NodeSimple.h"
#include "nuto/mechanics/elements/ElementInterface.h"
#include "nuto/mechanics/iga/Nurbs.h"
#include "nuto/mechanics/cell/Matrix.h"
#include "nuto/math/shapes/Spline.h"

namespace NuTo
{
template <int TDimParameter>
class ElementIga : public ElementInterface
{
public:
    ElementIga(const std::array<int, TDimParameter>& knotIDs, const Nurbs<TDimParameter>& NurbsGeometry)
        : mKnotIDs(knotIDs)
        , mNurbsGeometry(NurbsGeometry)
    {}

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
    //! @remark virtual to make it testable
    virtual NodeValues ExtractNodeValues(int instance = 0) const override
    {
        (void)instance;
        return NurbsGeometry().GetControlPointCoordinatesElement(mKnotIDs);
    }

    NMatrix GetNMatrix(NaturalCoords ipCoords) const override
    {
        return NuTo::Matrix::N(GetShapeFunctions(ipCoords), GetNumNodes(), GetDofDimension());
    }

    ShapeFunctions GetShapeFunctions(NaturalCoords ipCoords) const override
    {
        return NurbsGeometry().BasisFunctionsAndDerivativesRational(0, Transformation(ipCoords));
    }

    DerivativeShapeFunctionsNatural GetDerivativeShapeFunctions(NaturalCoords ipCoords) const override
    {
        return NurbsGeometry().BasisFunctionsAndDerivativesRational(1, Transformation(ipCoords));
    }

    int GetDofDimension() const override
    {
        return NurbsGeometry().GetDimension();
    }

    Eigen::VectorXi GetDofNumbering() const override
    {
        Eigen::VectorXi dofNumbering(GetNumNodes() * GetDofDimension());
        int i = 0;
        for (int iNode = 0; iNode < GetNumNodes(); ++iNode)
        {
            const auto& node = *(NurbsGeometry().GetControlPointElement(mKnotIDs, iNode));
            for (int iDof = 0; iDof < GetDofDimension(); ++iDof)
            {
                dofNumbering[i++] = node.GetDofNumber(iDof);
            }
        }
        return dofNumbering;
    }


    int GetNumNodes() const override
    {
        return NurbsGeometry().GetNumControlPointsElement();
    }

    const Shape& GetShape() const
    {
        return mShape;
    }

private:
    const Nurbs<TDimParameter>& NurbsGeometry() const
    {
        return mNurbsGeometry;
    }

    std::array<int, TDimParameter> mKnotIDs;
    std::reference_wrapper<const Nurbs<TDimParameter>> mNurbsGeometry;
    Spline mShape;

};
} /* NuTo */
