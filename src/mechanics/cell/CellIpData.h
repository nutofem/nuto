#pragma once

#include "mechanics/dofs/DofContainer.h"
#include "mechanics/cell/CellData.h"
#include "mechanics/cell/Jacobian.h"
#include "mechanics/cell/GradientOperators.h"

namespace NuTo
{

//! named pair for cellId and ipId to make the argument list shorter and avoid accidental mixup of both
struct Ids
{
    int cellId;
    int ipId;
};

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

    Ids GetIpId() const
    {
        return {mCellData.GetCellId(), mIpId};
    }

    Eigen::VectorXd Value(DofType dofType) const
    {
        return N(dofType) * mCellData.GetNodeValues(dofType);
    }

    Eigen::VectorXd Gradient(DofType dofType, const B::Interface& b = B::Gradient())
    {
        return B(dofType, b) * mCellData.GetNodeValues(dofType);
    }

    NMatrix N(const DofType& dofType) const
    {
        return mCellData.Elements().DofElement(dofType).GetNMatrix(mIPCoords);
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
