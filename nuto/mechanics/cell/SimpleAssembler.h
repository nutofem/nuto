#pragma once

#include "nuto/base/Group.h"
#include "nuto/mechanics/cell/CellInterface.h"
#include "nuto/mechanics/dofs/DofNumbering.h"
#include "nuto/mechanics/dofs/DofVector.h"

namespace NuTo
{
class SimpleAssembler
{
public:
    SimpleAssembler() = default;
    SimpleAssembler(DofInfo dofInfo);

    DofVector<double> BuildVector(const Group<CellInterface>& cells, std::vector<DofType> dofTypes,
                                CellInterface::VectorFunction f) const;

    Eigen::SparseMatrix<double> BuildMatrix(const Group<CellInterface>& cells, std::vector<DofType> dofTypes,
                                                         DofInfo dofInfo,
                                      CellInterface::MatrixFunction f) const;

    //! @brief Assembles a diagonally lumped matrix from local matrices calculated by f
    //! @param cells group of cells to be used for assembly
    //! @param dofTypes vector of dofTypes
    //! @return a dof vector that represents a collection of diagonal matrices
    //! @remark HRZ lumping is used here, the assumptions made are:
    //! - shape functions sum to 1, only then the total mass for a cell can be calculated by
    //!   summing over all entries of the local mass matrix
    //! - the total mass of a cell is the same for all components of the considered dof
    DofVector<double> BuildDiagonallyLumpedMatrix(const Group<CellInterface>& cells, std::vector<DofType> dofTypes,
                                                CellInterface::MatrixFunction f) const;

    void SetDofInfo(DofInfo dofInfo);

private:
    DofVector<double> ProperlyResizedVector(std::vector<DofType> dofTypes) const;

    DofInfo mDofInfo;

    void ThrowOnZeroDofNumbering(std::vector<DofType> dofTypes) const;
};
} /* NuTo */
