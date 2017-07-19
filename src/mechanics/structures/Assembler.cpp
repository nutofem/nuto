#include "mechanics/structures/Assembler.h"
#include "mechanics/nodes/NodeEnum.h"
#include "base/Exception.h"
#include "mechanics/nodes/NodeBase.h"
#include "math/SparseMatrixCSRVector2General.h"
#include "mechanics/structures/StructureOutputBlockVector.h"
#include "mechanics/dofSubMatrixStorage/BlockFullVector.h"

NuTo::Assembler::Assembler()
    : mNodeVectorChanged(false)
    , mConstraintMatrix(mDofStatus, false)
    , mConstraintMappingRhs(mDofStatus, false)
    , mConstraintRhs(mDofStatus)
{
}

void NuTo::Assembler::BuildGlobalDofs(const std::vector<NodeBase*>& rNodes)
{
    std::map<Node::eDof, int> numDofsMap;

    // build initial node numbering

    //
    // number Lagrange multipliers in constraint equations defined in StructureBase
    // currently removed
    // ConstraintNumberGlobalDofs(this->mNumDofs);
    //

    for (auto node : rNodes)
    {
        for (auto dof : node->GetDofTypes())
        {
            if (not node->IsDof(dof))
                continue;
            for (int i = 0; i < node->GetNum(dof); ++i)
            {
                node->SetDofNumber(dof, i, numDofsMap[dof]++);
            }
        }
    }

    mConstraintMatrix.AllocateSubmatrices();
    mConstraintMappingRhs.AllocateSubmatrices();
    mConstraintRhs.AllocateSubvectors();

    for (auto dof : mDofStatus.GetDofTypes())
    {
        int numConstraints = GetConstraints().GetNumEquations(dof);
        mDofStatus.SetNumDependentDofs(dof, numConstraints);
        mDofStatus.SetNumActiveDofs(dof, numDofsMap[dof] - numConstraints);
    }

    for (auto dof : mDofStatus.GetDofTypes())
    {
        auto& constraintMatrix = mConstraintMatrix(dof, dof);
        constraintMatrix = GetConstraints().BuildConstraintMatrix(dof, mDofStatus.GetNumDofs(dof));

        const int numActiveDofs = mDofStatus.GetNumActiveDofs(dof);
        const int numDependentDofs = mDofStatus.GetNumDependentDofs(dof);
        const int numDofs = numActiveDofs + numDependentDofs;

        // init RhsMatrix as a diagonal identity matrix
        auto& constraintMappingRhs = mConstraintMappingRhs(dof, dof);

        constraintMappingRhs.Resize(numDependentDofs, numDependentDofs);
        for (int i = 0; i < numDependentDofs; ++i)
            constraintMappingRhs.AddValue(i, i, 1.);

        // perform gauss algorithm
        std::vector<int> mappingInitialToNewOrdering;
        std::vector<int> mappingNewToInitialOrdering;

        constraintMatrix.Gauss(constraintMappingRhs, mappingNewToInitialOrdering, mappingInitialToNewOrdering);

        // move dependent dofs at the end
        // Warning!!! after this loop mappingNewToInitialOrdering is no longer valid !!!
        std::vector<int> tmpMapping;
        for (int dependentDofCount = 0; dependentDofCount < numDependentDofs; dependentDofCount++)
        {
            tmpMapping.push_back(numActiveDofs + dependentDofCount);
            mappingInitialToNewOrdering[mappingNewToInitialOrdering[dependentDofCount]] += numActiveDofs;
        }
        for (int activeDofCount = numDependentDofs; activeDofCount < numDofs; activeDofCount++)
        {
            tmpMapping.push_back(activeDofCount - numDependentDofs);
            mappingInitialToNewOrdering[mappingNewToInitialOrdering[activeDofCount]] -= numDependentDofs;
        }
        mappingNewToInitialOrdering.clear();

        // reorder columns
        constraintMatrix.ReorderColumns(tmpMapping);


        // remove columns of dependent dofs
        // check if the submatrix which is removed is a diagonal matrix
        const auto& columns = constraintMatrix.GetColumns();

        for (unsigned int iRow = 0; iRow < columns.size(); iRow++)
            for (unsigned int iPos = 0; iPos < columns[iRow].size(); iPos++)
            {
                int column = columns[iRow][iPos];
                if (column > numActiveDofs)
                    if (column - numActiveDofs != (int)iRow)
                        throw Exception(__PRETTY_FUNCTION__, "invalid matrix structure.");
            }


        constraintMatrix.RemoveLastColumns(numDependentDofs);

        // renumber dofs
        for (auto node : rNodes)
        {
            for (int i = 0; i < node->GetNum(dof); ++i)
            {
                int initialDofNumber = node->GetDof(dof, i);
                int newDofNumber = mappingInitialToNewOrdering[initialDofNumber];
                node->SetDofNumber(dof, i, newDofNumber);
            }
        }
    }

    // since only the diagonals were set, the off-diagonal submatrices have to be resized
    // to guarantee the right dimensions in arithmetic operations
    mConstraintMatrix.FixOffDiagonalDimensions();
    mConstraintMappingRhs.FixOffDiagonalDimensions();

    mConstraintMatrix.CheckDimensions();
    mConstraintMappingRhs.CheckDimensions();


    mNodeVectorChanged = false;
    GetConstraints().SetHaveChanged(false);

    // Build the Rhs once at the global time 0.
    // A call to GetRhsAfterGaussElimination would otherwise return an empty RHS.
    // This is required for static solutions that do not use a time integration
    // scheme. In those cases the RHS is not continuously updated and needs to
    // be calculated at some point. This point is here.
    ConstraintUpdateRhs(0);
}

NuTo::BlockFullVector<double> NuTo::Assembler::BuildRhsBeforeGaussElimination(double time) const
{
    ThrowIfRenumberingRequred();

    NuTo::BlockFullVector<double> rhsBeforeGaussElimination(mDofStatus);

    for (auto dof : mDofStatus.GetDofTypes())
        rhsBeforeGaussElimination[dof] = GetConstraints().GetRhs(dof, time);

    return rhsBeforeGaussElimination;
}

void NuTo::Assembler::ConstraintUpdateRhs(double time)
{
    ThrowIfRenumberingRequred();

    BlockFullVector<double> rhsBeforeGaussElimination = BuildRhsBeforeGaussElimination(time);

    // calculate the rhs vector of the constraint equations after the Gauss elimination using the mapping matrix
    mConstraintRhs = mConstraintMappingRhs * rhsBeforeGaussElimination;
}

void NuTo::Assembler::ThrowIfRenumberingRequred() const
{
    if (RenumberingRequired())
        throw Exception(__PRETTY_FUNCTION__, "build global numbering first");
}


NuTo::BlockFullVector<double> NuTo::Assembler::ApplyCMatrix(const StructureOutputBlockVector& vec, const BlockSparseMatrix& cMat)
{
    auto result = vec.J;

    if (not vec.J.GetDofStatus().HasInteractingConstraints())
        return result;

    for (auto dof : vec.J.GetDofStatus().GetActiveDofTypes())
        result[dof] -= cMat(dof, dof).TransMult(vec.K[dof]);

    return result;
}
