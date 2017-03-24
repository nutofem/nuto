#include "mechanics/structures/Assembler.h"
#include "mechanics/nodes/NodeEnum.h"
#include "mechanics/MechanicsException.h"
#include "mechanics/nodes/NodeBase.h"
#include "math/SparseMatrixCSRVector2General.h"
#include "mechanics/constraints/ConstraintLinear.h"

NuTo::Assembler::Assembler()
    : mNodeNumberingRequired(false)
    , mConstraintMatrix(mDofStatus, false)
    , mConstraintMappingRHS(mDofStatus, false)
    , mConstraintRHS(mDofStatus)
{
}

void NuTo::Assembler::BuildGlobalDofs(boost::ptr_map<int, NodeBase>& rNodes)
{
    mNodeNumberingRequired = false;
    std::map<Node::eDof, int> numDofsMap;

    // build initial node numbering

    //
    // number Lagrange multipliers in constraint equations defined in StructureBase
    // currently removed
    // ConstraintNumberGlobalDofs(this->mNumDofs);
    //

    for (auto it : rNodes)
    {
        NodeBase& node = *it->second;
        for(auto dof : node.GetDofTypes())
        {
            if (not node.IsDof(dof))
                continue;
            for (int i = 0; i < node.GetNum(dof); ++i)
            {
                node.SetDofNumber(dof, i, numDofsMap[dof]++);
            }
        }
    }

    mConstraintMatrix.AllocateSubmatrices();
    mConstraintMappingRHS.AllocateSubmatrices();
    mConstraintRHS.AllocateSubvectors();


    for (auto dof : mDofStatus.GetDofTypes())
    {
        int numConstraints = ConstraintGetNumLinearConstraints(dof);
        mDofStatus.SetNumDependentDofs(dof, numConstraints);
        mDofStatus.SetNumActiveDofs(dof, numDofsMap[dof] - numConstraints);
    }


    // build constraint matrix for all dofs
    mConstraintMatrix = ConstraintGetConstraintMatrixBeforeGaussElimination();

    for (auto dof : mDofStatus.GetDofTypes())
    {
        auto& constraintMatrix = mConstraintMatrix(dof, dof);

        const int numActiveDofs = mDofStatus.GetNumActiveDofs(dof);
        const int numDependentDofs = mDofStatus.GetNumDependentDofs(dof);
        const int numDofs = numActiveDofs + numDependentDofs;

        // init RHSMatrix as a diagonal identity matrix
        auto& constraintMappingRHS = mConstraintMappingRHS(dof, dof);

        constraintMappingRHS.Resize(numDependentDofs, numDependentDofs);
        for (int i = 0; i < numDependentDofs; ++i)
            constraintMappingRHS.AddValue(i, i, 1.);

        // perform gauss algorithm
        std::vector<int> mappingInitialToNewOrdering;
        std::vector<int> mappingNewToInitialOrdering;

        constraintMatrix.Gauss(constraintMappingRHS, mappingNewToInitialOrdering, mappingInitialToNewOrdering);

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

        for (auto it : rNodes)
        {
            NodeBase& node = *it->second;
            for (int i = 0; i < node.GetNum(dof); ++i)
            {
                int initialDofNumber = node.GetDof(dof, i);
                int newDofNumber = mappingInitialToNewOrdering[initialDofNumber];
                node.SetDofNumber(dof, i, newDofNumber);
            }
        }
    }

    // since only the diagonals were set, the off-diagonal submatrices have to be resized
    // to guarantee the right dimensions in arithmetic operations
    mConstraintMatrix.FixOffDiagonalDimensions();
    mConstraintMappingRHS.FixOffDiagonalDimensions();

    mConstraintMatrix.CheckDimensions();
    mConstraintMappingRHS.CheckDimensions();


    // number Lagrange multipliers in constraint equations defined in StructureBase
    // currently removed
    // renumber DOFS in constraints (Lagrange multiplier)
    //        ConstraintRenumberGlobalDofs(mappingInitialToNewOrdering);

    // calculate current rhs matrix
    ConstraintUpdateRHSAfterGaussElimination();
    mNodeNumberingRequired = false;
}


int NuTo::Assembler::ConstraintGetNumLinearConstraints(std::string rDof) const

{
    return ConstraintGetNumLinearConstraints(Node::DofToEnum(rDof));
}

int NuTo::Assembler::ConstraintGetNumLinearConstraints(Node::eDof rDof) const
{
    int numLinearConstraints = 0;
    for (auto itConstraint : mConstraintMap)
    {
        const auto& constraint = itConstraint.second;
        if (constraint->GetDofType() == rDof)
            numLinearConstraints += constraint->GetNumLinearConstraints();
    }
    return numLinearConstraints;
}

NuTo::BlockSparseMatrix NuTo::Assembler::ConstraintGetConstraintMatrixBeforeGaussElimination() const
{
    if (mNodeNumberingRequired)
    {
        throw MechanicsException(__PRETTY_FUNCTION__, "build global numbering first");
    }

    BlockSparseMatrix constraintMatrix(mDofStatus, false);

    for (auto dofType : mDofStatus.GetDofTypes())
    {
        int numLinearConstraints = ConstraintGetNumLinearConstraints(dofType);
        int curConstraintEquations = 0;

        constraintMatrix(dofType, dofType).Resize(numLinearConstraints, mDofStatus.GetNumDofs(dofType));

        for (auto itConstraint : mConstraintMap)
        {
            if (itConstraint->second->GetNumLinearConstraints() > 0)
            {
                try
                {
                    if (itConstraint->second->GetDofType() == dofType)
                        itConstraint->second->AsConstraintLinear()->AddToConstraintMatrix(
                                curConstraintEquations, constraintMatrix(dofType, dofType));
                }
                catch (MechanicsException& e)
                {
                    e.AddMessage(__PRETTY_FUNCTION__, "mechanics exception while building constraint matrix for "
                                                      "constraint with nonzero number of linear components.");
                    throw;
                }
                catch (...)
                {
                    throw MechanicsException(__PRETTY_FUNCTION__, "error building constraint matrix for constraint "
                                                                  "with nonzero number of linear components.");
                }
            }
        }

        if (curConstraintEquations != numLinearConstraints)
        {
            std::cout << "curConstraintEquations " << curConstraintEquations << std::endl;
            std::cout << "numConstraintEquations " << numLinearConstraints << std::endl;
            throw MechanicsException(__PRETTY_FUNCTION__,
                                     "Internal error, there is something wrong with the constraint equations.");
        }
    }
    return constraintMatrix;
}

NuTo::BlockFullVector<double> NuTo::Assembler::ConstraintGetRHSBeforeGaussElimination()
{
    if (mNodeNumberingRequired)
    {
        throw MechanicsException(__PRETTY_FUNCTION__, "build global numbering first");
    }

    NuTo::BlockFullVector<double> rhsBeforeGaussElimination(mDofStatus);

    for (auto dof : mDofStatus.GetDofTypes())
    {

        int numLinearConstraints = ConstraintGetNumLinearConstraints(dof);

        rhsBeforeGaussElimination[dof].resize(numLinearConstraints);

        // calculate the rhs vector of the constraint equations before the Gauss elimination
        int curConstraintEquations = 0;
        for (auto itConstraint : mConstraintMap)
        {
            if (itConstraint.second->GetNumLinearConstraints() > 0)
            {
                try
                {
                    if (itConstraint.second->GetDofType() == dof)
                        itConstraint.second->AsConstraintLinear()->GetRHS(curConstraintEquations,
                                                                          rhsBeforeGaussElimination[dof]);
                }
                catch (MechanicsException& e)
                {
                    e.AddMessage(__PRETTY_FUNCTION__,
                                 "mechanics exception while building rhs vector after gauss elimination.");
                    throw;
                }
                catch (...)
                {
                    throw MechanicsException(__PRETTY_FUNCTION__,
                                             "mechanics exception while building rhs vector after gauss elimination.");
                }
            }
        }

        if (curConstraintEquations != numLinearConstraints)
        {
            std::cout << "curConstraintEquations " << curConstraintEquations << std::endl;
            std::cout << "numConstraintEquations " << numLinearConstraints << std::endl;
            throw MechanicsException(__PRETTY_FUNCTION__,
                                     "Internal error, there is something wrong with the constraint equations.");
        }
    }
    return rhsBeforeGaussElimination;
}

void NuTo::Assembler::ConstraintUpdateRHSAfterGaussElimination()
{
    if (mNodeNumberingRequired)
    {
        throw MechanicsException(__PRETTY_FUNCTION__, "build global numbering first");
    }

    BlockFullVector<double> rhsBeforeGaussElimination = ConstraintGetRHSBeforeGaussElimination();

    if (mConstraintMappingRHS.GetNumColumns() != rhsBeforeGaussElimination.GetNumRows())
    {
        throw MechanicsException(__PRETTY_FUNCTION__, "here is something wrong in the implementation.");
    }

    // calculate the rhs vector of the constraint equations after the Gauss elimination using the mapping matrix
    mConstraintRHS = mConstraintMappingRHS * rhsBeforeGaussElimination;
}

double NuTo::Assembler::ConstraintGetRHS(int rConstraintEquation) const
{
    auto it = mConstraintMap.find(rConstraintEquation);
    if (it == mConstraintMap.end())
    {
        throw MechanicsException(__PRETTY_FUNCTION__, "Constraint equation does not exist.");
    }
    return it->second->GetRHS();
}
