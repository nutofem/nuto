#pragma once

#include "base/Group.h"
#include "mechanics/cell/CellInterface.h"
#include "mechanics/dofs/DofNumbering.h"
#include "mechanics/dofs/GlobalDofVector.h"
#include "mechanics/dofs/GlobalDofMatrixSparse.h"

namespace NuTo
{
class SimpleAssembler
{
public:
    SimpleAssembler() = default;
    SimpleAssembler(DofInfo dofInfo);

    GlobalDofVector BuildVector(const Group<CellInterface>& cells, std::vector<DofType> dofTypes,
                                CellInterface::VectorFunction f) const;

    GlobalDofMatrixSparse BuildMatrix(const Group<CellInterface>& cells, std::vector<DofType> dofTypes,
                                      CellInterface::MatrixFunction f) const;

    //! @brief Assembles a diagonally lumped matrix from local matrices calculated by f
    //! @param cells group of cells to be used for assembly
    //! @param dofTypes vector of dofTypes
    //! @return a global dof vector that represents a collection of diagonal matrices
    //! @remark HRZ lumping is used here, the assumptions made are:
    //! - shape functions sum to 1, only then the total mass for a cell can be calculated by
    //!   summing over all entries of the local mass matrix
    //! - the total mass of a cell is the same for all components of the considered dof
    GlobalDofVector BuildDiagonallyLumpedMatrix(const Group<CellInterface>& cells, std::vector<DofType> dofTypes,
                                                CellInterface::MatrixFunction f) const;

    void SetDofInfo(DofInfo dofInfo);

private:
    GlobalDofVector ProperlyResizedGlobalVector(std::vector<DofType> dofTypes) const;
    GlobalDofMatrixSparse ProperlyResizedGlobalMatrix(std::vector<DofType> dofTypes) const;

    DofInfo mDofInfo;

    void ThrowOnZeroDofNumbering(std::vector<DofType> dofTypes) const;
};
} /* NuTo */
