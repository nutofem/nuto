#pragma once

#include "mechanics/dofs/DofContainer.h"
#include "mechanics/cell/CellData.h"
#include "mechanics/cell/Jacobian.h"
#include "mechanics/cell/GradientOperators.h"
#include "mechanics/cell/CellIds.h"

namespace NuTo
{

//! @brief Similar to NuTo::CellData, but for N and B
class CellIpData
{
public:
    CellIpData(const CellData& cellData, NuTo::Jacobian jacobian, NaturalCoords ipCoords, int ipId)
        : mCellData(cellData)
        , mJacobian(std::move(jacobian))
        , mIPCoords(std::move(ipCoords))
        , mIpId(ipId)
    {
    }

    Eigen::VectorXd GlobalCoordinates() const
    {
        return Interpolate(mCellData.Elements().CoordinateElement(), mIPCoords);
    }

    CellIds Ids() const
    {
        return {mCellData.GetCellId(), mIpId};
    }

    Eigen::VectorXd Value(DofType dofType) const
    {
        return N(dofType) * mCellData.GetNodeValues(dofType);
    }

    Eigen::VectorXd Gradient(DofType dofType, const B::Interface& b = B::Gradient()) const
    {
        return B(dofType, b) * mCellData.GetNodeValues(dofType);
    }

    const Eigen::MatrixXd& N(const DofType& dofType) const
    {
        Eigen::MatrixXd& N = mNs[dofType];
        if (N.size() == 0) // simplest memoization using a mutable mNs to keep it const
            N = mCellData.Elements().DofElement(dofType).GetNMatrix(mIPCoords);
        return N;
    }

    const Eigen::MatrixXd& B(DofType dofType, const B::Interface& b = B::Gradient()) const
    {
        Eigen::MatrixXd& B = mBs[dofType];
        if (B.size() == 0) // simplest memoization using a mutable mBs to keep it const
            B = b(CalculateDerivativeShapeFunctionsGlobal(dofType));
        return B;
    }

private:
    DerivativeShapeFunctionsGlobal CalculateDerivativeShapeFunctionsGlobal(const DofType& dofType) const
    {
        DerivativeShapeFunctionsNatural dShapeNatural =
                mCellData.Elements().DofElement(dofType).GetDerivativeShapeFunctions(mIPCoords);
        return mJacobian.TransformDerivativeShapeFunctions(dShapeNatural);
    }

    const CellData& mCellData;
    NuTo::Jacobian mJacobian;
    NaturalCoords mIPCoords;
    int mIpId;

    mutable DofContainer<Eigen::MatrixXd> mBs;
    mutable DofContainer<Eigen::MatrixXd> mNs;
};
} /* NuTo */
