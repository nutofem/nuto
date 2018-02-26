#pragma once

#include "mechanics/dofs/DofContainer.h"
#include "mechanics/cell/CellData.h"
#include "mechanics/cell/Jacobian.h"
#include "mechanics/cell/DifferentialOperators.h"
#include "mechanics/cell/CellIds.h"

namespace NuTo
{

//! @brief Similar to NuTo::CellData. Memoizes N and B and allows access to the values and gradients of a dof at this
//! integration point
class CellIpData
{
public:
    //! ctor
    //! @param cellData cell data that is constant for all integration points of a cell (cellId, node values)
    //! @param jacobian transformation from the natural (xi, eta,...) system to the global system (x,y,..)
    //! @param ipCoords coordinates of the current integration point
    //! @param ipId id of the current integration point
    CellIpData(const CellData& cellData, NuTo::Jacobian jacobian, NaturalCoords ipCoords, int ipId)
        : mCellData(cellData)
        , mJacobian(std::move(jacobian))
        , mIPCoords(std::move(ipCoords))
        , mIpId(ipId)
    {
    }

    //! Caluclate the global integration point coordinates
    //! @return global integration point coordinates
    Eigen::VectorXd GlobalCoordinates() const
    {
        return Interpolate(mCellData.Elements().CoordinateElement(), mIPCoords);
    }

    //! Access to the cellId and ipId, compressed in CellIds
    //! @return named pair of cellId and ipId
    CellIds Ids() const
    {
        return {mCellData.GetCellId(), mIpId};
    }

    //! Calculates the value of a dof at the integration point
    //! @param dofType dof type that is interpolated
    //! @return interpolated dof value
    Eigen::VectorXd Value(DofType dofType) const
    {
        return N(dofType) * NodeValueVector(dofType);
    }
  
    const NuTo::Jacobian& GetJacobian() const
    {
        return mJacobian;
    }

    double Value(ScalarDofType dofType) const
    {
        return Value(DofType(dofType))[0];
    }

    //! Calculates the gradient (derivative of the value with respect to x) for a given dof type at the integration
    //! point
    //! @param dofType dof type that is evaluated
    //! @param b gradient operator that determines how to calculate the derivative (e.g. B::Gradient() for scalars or
    //! B::Strain() for engineering strains)
    //! @return gradient of the dof
    Eigen::VectorXd Apply(DofType dofType, const Nabla::Interface& b) const
    {
        return B(dofType, b) * NodeValueVector(dofType);
    }

    //! Returns a memoized copy of the N matrix for a given dof type
    //! @param dofType
    //! @return N matrix
    const Eigen::MatrixXd& N(DofType dofType) const
    {
        Eigen::MatrixXd& N = mNs[dofType];
        if (N.size() == 0) // simplest memoization using a mutable mNs to keep it const
            N = mCellData.Elements().DofElement(dofType).GetNMatrix(mIPCoords);
        return N;
    }

    //! Returns a memoized copy of the B matrix for a given dof type
    //! @param dofType
    //! @param b gradient operator that determines how to calculate the derivative (e.g. B::Gradient() for scalars or
    //! B::Strain() for engineering strains)
    //! @return B matrix
    const Eigen::MatrixXd& B(DofType dofType, const Nabla::Interface& b) const
    {
        Eigen::MatrixXd& B = mBs[dofType];
        if (B.size() == 0) // simplest memoization using a mutable mBs to keep it const
            B = b(CalculateDerivativeShapeFunctionsGlobal(dofType));
        return B;
    }

    //! Returns memoized nodal values
    //! @param dofType
    //! @return nodal values for dofType
    const NodeValues& NodeValueVector(DofType dofType) const
    {
        return mCellData.GetNodeValues(dofType);
    }

private:
    //! Transforms the derivative shape functions from the natural coordinate system (dN_d(xi, eta, ...)) to the global
    //! coordinate system (dN_d(x,y,...))
    DerivativeShapeFunctionsGlobal CalculateDerivativeShapeFunctionsGlobal(DofType dofType) const
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
