#include "mechanics/structures/Assembler.h"
#include "mechanics/nodes/NodeEnum.h"
#include "mechanics/MechanicsException.h"
#include "mechanics/nodes/NodeBase.h"
#include "math/SparseMatrixCSRVector2General.h"
#include "mechanics/constraints/ConstraintLinear.h"

NuTo::Assembler::Assembler()
    : mNodeNumberingRequired(false)
    , mConstraintMatrix(mDofStatus, false)
    , mConstraintMappingRhs(mDofStatus, false)
    , mConstraintRhs(mDofStatus)
{
}

void NuTo::Assembler::BuildGlobalDofs(const std::vector<NodeBase*>& rNodes)
{
    mNodeNumberingRequired = false;
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
        int numConstraints = mConstraints.GetNumEquations(dof);
        mDofStatus.SetNumDependentDofs(dof, numConstraints);
        mDofStatus.SetNumActiveDofs(dof, numDofsMap[dof] - numConstraints);
    }

    std::cout << mDofStatus << std::endl;

    for (auto dof : mDofStatus.GetDofTypes())
    {
        auto& constraintMatrix = mConstraintMatrix(dof, dof);
        constraintMatrix.Resize(mConstraints.GetNumEquations(dof), mDofStatus.GetNumDofs(dof));
        mConstraints.BuildConstraintMatrix(constraintMatrix, dof);

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
                        throw MechanicsException(__PRETTY_FUNCTION__, "invalid matrix structure.");
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


    // number Lagrange multipliers in constraint equations defined in StructureBase
    // currently removed
    // renumber DOFS in constraints (Lagrange multiplier)
    //        ConstraintRenumberGlobalDofs(mappingInitialToNewOrdering);

    // calculate current rhs matrix
    mNodeNumberingRequired = false;
}


int NuTo::Assembler::ConstraintGetNumLinearConstraints(std::string rDof) const
{
    return ConstraintGetNumLinearConstraints(Node::DofToEnum(rDof));
}

int NuTo::Assembler::ConstraintGetNumLinearConstraints(Node::eDof rDof) const
{
    return mConstraints.GetNumEquations(rDof);
}


NuTo::BlockFullVector<double> NuTo::Assembler::ConstraintGetRhsBeforeGaussElimination(double time) const
{
    if (mNodeNumberingRequired)
        throw MechanicsException(__PRETTY_FUNCTION__, "build global numbering first");

    NuTo::BlockFullVector<double> rhsBeforeGaussElimination(mDofStatus);

    for (auto dof : mDofStatus.GetDofTypes())
        rhsBeforeGaussElimination[dof] = mConstraints.GetRhs(dof, time);
    
    return rhsBeforeGaussElimination;
}

void NuTo::Assembler::ConstraintUpdateRhs(double time)
{
    if (mNodeNumberingRequired)
        throw MechanicsException(__PRETTY_FUNCTION__, "build global numbering first");

    BlockFullVector<double> rhsBeforeGaussElimination = ConstraintGetRhsBeforeGaussElimination(time);

    // calculate the rhs vector of the constraint equations after the Gauss elimination using the mapping matrix
    mConstraintRhs = mConstraintMappingRhs * rhsBeforeGaussElimination;
}

void NuTo::Assembler::AddEquation(NuTo::Node::eDof dof, Constraint::Equation equation)
{
    mConstraints.AddEquation(dof, equation);
    mNodeNumberingRequired = true;
}

void NuTo::Assembler::AddEquations(NuTo::Node::eDof dof, std::vector<Constraint::Equation> equations)
{
    for (auto& equation : equations)
        AddEquation(dof, equation);
}
