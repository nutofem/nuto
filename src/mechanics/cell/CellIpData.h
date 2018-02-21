#pragma once

#include "mechanics/dofs/DofContainer.h"
#include "mechanics/elements/ElementCollection.h"
#include "mechanics/cell/Jacobian.h"
#include "mechanics/cell/GradientOperators.h"

namespace NuTo
{

//! @brief Similar to NuTo::CellData, but for N and B
class CellIpData
{
public:
    CellIpData(const ElementCollection& elements, NuTo::Jacobian jacobian, NaturalCoords ipCoords, int ipId)
        : mElements(elements)
        , mJacobian(std::move(jacobian))
        , mIPCoords(std::move(ipCoords))
        , mIpId(ipId)
    {
    }

    Eigen::VectorXd GlobalCoordinates() const
    {
        return Interpolate(mElements.CoordinateElement(), mIPCoords);
    }

    int GetIpId() const
    {
        return mIpId;
    }

    NMatrix N(const DofType& dofType) const
    {
        return mElements.DofElement(dofType).GetNMatrix(mIPCoords);
    }

    const Eigen::MatrixXd& B(DofType dofType, const B::Interface& b) const
    {
        if (not mBs.Has(dofType)) // simplest memoization using a mutable mBs to keep it const
            mBs[dofType] = b(CalculateDerivativeShapeFunctionsGlobal(dofType));

        return mBs[dofType];
    }

private:
    DerivativeShapeFunctionsGlobal CalculateDerivativeShapeFunctionsGlobal(const DofType& dofType) const
    {
        DerivativeShapeFunctionsNatural dShapeNatural =
                mElements.DofElement(dofType).GetDerivativeShapeFunctions(mIPCoords);
        return mJacobian.TransformDerivativeShapeFunctions(dShapeNatural);
    }

    const ElementCollection& mElements;
    NuTo::Jacobian mJacobian;
    NaturalCoords mIPCoords;
    int mIpId;

    mutable DofContainer<Eigen::MatrixXd> mBs;
};
} /* NuTo */
