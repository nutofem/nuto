#pragma once

#include "nuto/mechanics/dofs/DofContainer.h"
#include "nuto/mechanics/cell/CellData.h"
#include "nuto/mechanics/cell/Jacobian.h"
#include "nuto/mechanics/cell/DifferentialOperators.h"
#include "nuto/mechanics/cell/CellIds.h"

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
    //! @param instance id of the dof instance
    //! @return interpolated dof value
    Eigen::VectorXd Value(DofType dofType, int instance = 0) const
    {
        return N(dofType) * NodeValueVector(dofType, instance);
    }

    const NuTo::Jacobian& GetJacobian() const
    {
        return mJacobian;
    }

    double Value(ScalarDofType dofType, int instance = 0) const
    {
        return Value(DofType(dofType), instance)[0];
    }

    //! Calculates the gradient (derivative of the value with respect to x) for a given dof type at the integration
    //! point
    //! @param dofType dof type that is evaluated
    //! @param b gradient operator that determines how to calculate the derivative (e.g. B::Gradient() for scalars or
    //! @param instance id of the dof instance
    //! B::Strain() for engineering strains)
    //! @return gradient of the dof
    Eigen::VectorXd Apply(DofType dofType, const Nabla::Interface& b, int instance = 0) const
    {
        return B(dofType, b) * NodeValueVector(dofType, instance);
    }

    //! Returns the N matrix for a given dof type
    //! @return N matrix
    const Eigen::MatrixXd N(DofType dofType) const
    {
        return mCellData.Elements().DofElement(dofType).GetNMatrix(mIPCoords);
    }

    //! Returns the B matrix for a given dof type
    //! @param b gradient operator that determines how to calculate the derivative (e.g. B::Gradient() for scalars or
    //! B::Strain() for engineering strains)
    //! @return B matrix
    const Eigen::MatrixXd B(DofType dofType, const Nabla::Interface& b) const
    {
        return b(CalculateDerivativeShapeFunctionsGlobal(dofType));
    }

    //! Returns memoized nodal values
    //! @return nodal values for dofType
    //! @param instance id of the dof instance
    const Eigen::VectorXd& NodeValueVector(DofType dofType, int instance = 0) const
    {
        return mCellData.GetNodeValues(dofType, instance);
    }

private:
    //! Transforms the derivative shape functions from the natural coordinate system (dN_d(xi, eta, ...)) to the global
    //! coordinate system (dN_d(x,y,...))
    Eigen::MatrixXd CalculateDerivativeShapeFunctionsGlobal(DofType dofType) const
    {
        Eigen::MatrixXd dShapeNatural = mCellData.Elements().DofElement(dofType).GetDerivativeShapeFunctions(mIPCoords);
        return mJacobian.TransformDerivativeShapeFunctions(dShapeNatural);
    }

    const CellData& mCellData;
    NuTo::Jacobian mJacobian;
    NaturalCoords mIPCoords;
    int mIpId;
};
} /* NuTo */
