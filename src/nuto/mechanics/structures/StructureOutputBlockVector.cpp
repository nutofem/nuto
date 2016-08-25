
#include "nuto/math/FullVector.h"
#include "nuto/math/SparseMatrixCSRVector2.h"
#include "nuto/mechanics/structures/StructureOutputBlockVector.h"
#include "nuto/mechanics/dofSubMatrixStorage/BlockSparseMatrix.h"
#include "nuto/mechanics/dofSubMatrixStorage/DofStatus.h"

NuTo::StructureOutputBlockVector::StructureOutputBlockVector(const DofStatus& rDofStatus, bool rAutomaticResize)
    : StructureOutputBase(),
      J(rDofStatus),
      K(rDofStatus)
{
    if (rAutomaticResize)
        Resize(rDofStatus.GetNumActiveDofsMap(), rDofStatus.GetNumDependentDofsMap());
}

NuTo::StructureOutputBlockVector::~StructureOutputBlockVector()
{

}

NuTo::StructureOutputBlockVector &NuTo::StructureOutputBlockVector::operator=(const NuTo::StructureOutputBlockVector &rOther)
{
    J = rOther.J;
    K = rOther.K;
    return *this;
}


void NuTo::StructureOutputBlockVector::AddElementVector(
        const NuTo::BlockFullVector<double>& rElementVector,
        const NuTo::BlockFullVector<int>& rGlobalRowDofNumbers)
{
    const auto& activeDofTypes = J.GetDofStatus().GetActiveDofTypes();
    const auto& numActiveDofTypeMap = J.GetDofStatus().GetNumActiveDofsMap();
    for (auto dofRow : activeDofTypes)
    {
        const auto& elementVector = rElementVector[dofRow];
        const auto& globalRowDofs = rGlobalRowDofNumbers[dofRow];

        int numActiveDofsRow = numActiveDofTypeMap.at(dofRow);

        assert(elementVector.GetNumRows() == globalRowDofs.GetNumRows());

        auto& activeRow = J[dofRow];
        auto& dependentRow = K[dofRow];

        for (int iRow = 0; iRow < globalRowDofs.GetNumRows(); ++iRow)
        {
            int globalRowDof = globalRowDofs[iRow];
            if (globalRowDof < numActiveDofsRow)
            {
                activeRow[globalRowDof] += elementVector[iRow];
            }
            else
            {
                dependentRow[globalRowDof - numActiveDofsRow] += elementVector[iRow];
            }
        }
    }

}

void NuTo::StructureOutputBlockVector::ApplyCMatrix(BlockFullVector<double>& rResidual, const BlockSparseMatrix& rCmat) const
{
    rResidual = J;

    if (not J.GetDofStatus().HasInteractingConstraints())
        return;

    for (auto dof : J.GetDofStatus().GetActiveDofTypes())
        rResidual[dof] -= rCmat(dof,dof).TransMult(K[dof]);
}

void NuTo::StructureOutputBlockVector::ApplyCMatrix(const BlockSparseMatrix& rCmat)
{
    if (not J.GetDofStatus().HasInteractingConstraints())
        return;

    for (auto dof : J.GetDofStatus().GetActiveDofTypes())
        J[dof] += rCmat(dof,dof).TransMult(K[dof]);
}

void NuTo::StructureOutputBlockVector::Resize(const std::map<Node::eDof, int>& rNumActiveDofsMap, const std::map<Node::eDof, int>& rNumDependentDofsMap)
{
    assert(rNumActiveDofsMap.size() == rNumDependentDofsMap.size());

    J.Resize(rNumActiveDofsMap   );
    K.Resize(rNumDependentDofsMap);
}

NuTo::StructureOutputBlockVector& NuTo::StructureOutputBlockVector::operator +=(const StructureOutputBlockVector& rRhs)
{
    J += rRhs.J;
    K += rRhs.K;
    return *this;
}

NuTo::StructureOutputBlockVector& NuTo::StructureOutputBlockVector::operator -=(const StructureOutputBlockVector& rRhs)
{
    J -= rRhs.J;
    K -= rRhs.K;
    return *this;
}
NuTo::StructureOutputBlockVector& NuTo::StructureOutputBlockVector::operator *=(double rRhs)
{
    J *= rRhs;
    K *= rRhs;
    return *this;
}

NuTo::StructureOutputBlockVector& NuTo::StructureOutputBlockVector::operator /=(double rRhs)
{
    J /= rRhs;
    K /= rRhs;
    return *this;
}


//! @brief stream operator for outputs with cout or files
std::ostream& NuTo::operator<<(std::ostream &rOut, const NuTo::StructureOutputBlockVector& rStructureOutputBlockVector)
{
    rOut << "Active Dofs" << std::endl;
    rOut << rStructureOutputBlockVector.J << std::endl;
    rOut << "Dependent Dofs" << std::endl;
    rOut << rStructureOutputBlockVector.K << std::endl;
    return rOut;
}

void NuTo::StructureOutputBlockVector::SetZero()
{
    J.SetZero();
    K.SetZero();
}
