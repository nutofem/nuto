#include "nuto/mechanics/structures/StructureOutputBlockMatrix.h"


#include "nuto/math/SparseMatrixCSRVector2.h"
#include "nuto/mechanics/constitutive/ConstitutiveBase.h"
#include "nuto/mechanics/dofSubMatrixStorage/BlockFullMatrix.h"
#include "nuto/mechanics/dofSubMatrixStorage/BlockFullVector.h"
#include "nuto/mechanics/dofSubMatrixStorage/DofStatus.h"
#include "nuto/mechanics/structures/StructureOutputBlockVector.h"
#include "nuto/mechanics/elements/ElementBase.h"


NuTo::StructureOutputBlockMatrix::StructureOutputBlockMatrix(const DofStatus &rDofStatus, bool rAutomaticResize)
    : StructureOutputBase(),
      JJ(rDofStatus, true),
      JK(rDofStatus, false),
      KJ(rDofStatus, false),
      KK(rDofStatus, true)
{
    if (rAutomaticResize)
        Resize(rDofStatus.GetNumActiveDofsMap(), rDofStatus.GetNumDependentDofsMap());
}

NuTo::StructureOutputBlockMatrix::~StructureOutputBlockMatrix()
{}

void NuTo::StructureOutputBlockMatrix::AddElementMatrix(
        const ElementBase* rElementPtr,
        const NuTo::BlockFullMatrix<double>& rElementMatrix,
        const NuTo::BlockFullVector<int>& rGlobalRowDofNumbers,
        const NuTo::BlockFullVector<int>& rGlobalColumnDofNumbers,
        double rAddValueTolerance, bool rAssembleKJKK)
{
    const auto& activeDofTypes = JJ.GetDofStatus().GetActiveDofTypes();
    const auto& numActiveDofTypeMap = JJ.GetDofStatus().GetNumActiveDofsMap();

    for (auto dofRow : activeDofTypes)
        for (auto dofCol : activeDofTypes)
        {
            // TODO: Do not loop over all possible combinations of DOFs but over a list of combinations created by the constitutive law of the corresponding element. What if an element has multiple constitutive laws assigned?
            if (not rElementPtr->GetConstitutiveLaw(0).CheckDofCombinationComputable(dofRow, dofCol, 0))
                continue;

            const auto& elementMatrix = rElementMatrix(dofRow, dofCol);
            const auto& globalRowDofs = rGlobalRowDofNumbers[dofRow];
            const auto& globalColDofs = rGlobalColumnDofNumbers[dofCol];

            int numActiveDofsRow = numActiveDofTypeMap.at(dofRow);
            int numActiveDofsCol = numActiveDofTypeMap.at(dofCol);

//            std::cout << "elementMatrix.GetNumRows()    " << elementMatrix.GetNumRows()    << std::endl;
//            std::cout << "globalRowDofs.GetNumRows()    " << globalRowDofs.GetNumRows()    << std::endl;
//            std::cout << "elementMatrix.GetNumColumns() " << elementMatrix.GetNumColumns() << std::endl;
//            std::cout << "globalColDofs.GetNumRows()    " << globalColDofs.GetNumRows()    << std::endl;


            assert(elementMatrix.GetNumRows() == globalRowDofs.GetNumRows());
            assert(elementMatrix.GetNumColumns() == globalColDofs.GetNumRows());

            for (int iRow = 0; iRow < globalRowDofs.GetNumRows(); ++iRow)
            {
                int globalRowDof = globalRowDofs[iRow];
                if (globalRowDof < numActiveDofsRow)
                {
                    auto& activeCol = JJ(dofRow, dofCol);
                    auto& dependentCol = JK(dofRow, dofCol);

                    for (int iCol = 0; iCol < globalColDofs.GetNumRows(); ++iCol)
                    {
                        double value = elementMatrix(iRow, iCol);
                        if (std::abs(value - rAddValueTolerance) > 0. )
                        {
                            int globalColDof = globalColDofs[iCol];
                            if (globalColDof < numActiveDofsCol)
                                activeCol.AddValue(globalRowDof, globalColDof, value);
                            else
                                dependentCol.AddValue(globalRowDof, globalColDof - numActiveDofsCol, value);
                        }
                    }
                }
                else
                {
                    if (rAssembleKJKK)
                    {
                        auto& activeCol = KJ(dofRow, dofCol);
                        auto& dependentCol = KK(dofRow, dofCol);

                        globalRowDof -= numActiveDofsRow;

                        for (int iCol = 0; iCol < globalColDofs.GetNumRows(); ++iCol)
                        {
                            double value = elementMatrix(iRow, iCol);
                            if (std::abs(value - rAddValueTolerance) > 0. )
                            {
                                int globalColDof = globalColDofs[iCol];
                                if (globalColDof < numActiveDofsCol)
                                    activeCol.AddValue(globalRowDof, globalColDof, value);
                                else
                                    dependentCol.AddValue(globalRowDof, globalColDof - numActiveDofsCol, value);
                            }
                        }
                    }
                }
            }

        }
}

void NuTo::StructureOutputBlockMatrix::AddElementVectorDiagonal(
        const NuTo::BlockFullVector<double>& rElementVector,
        const NuTo::BlockFullVector<int>& rGlobalRowDofNumbers,
        double rAddValueTolerance)
{
    const auto& activeDofTypes = JJ.GetDofStatus().GetActiveDofTypes();
    const auto& numActiveDofTypeMap = JJ.GetDofStatus().GetNumActiveDofsMap();
    for (auto dofRow : activeDofTypes)
    {
        const auto& elementVector = rElementVector[dofRow];
        const auto& globalRowDofs = rGlobalRowDofNumbers[dofRow];

        int numActiveDofsRow = numActiveDofTypeMap.at(dofRow);

        assert(elementVector.GetNumRows() == globalRowDofs.GetNumRows());

        auto& activeRow = JJ(dofRow, dofRow);
        auto& dependentRow = KK(dofRow, dofRow);

        for (int iRow = 0; iRow < globalRowDofs.GetNumRows(); ++iRow)
        {
            double value = elementVector[iRow];
            if (std::abs(value - rAddValueTolerance) > 0. )
            {
                int globalRowDof = globalRowDofs[iRow];
                if (globalRowDof < numActiveDofsRow)
                    activeRow.AddValue(globalRowDof, globalRowDof, value);
                else
                {
                    globalRowDof -= numActiveDofsRow;
                    dependentRow.AddValue(globalRowDof, globalRowDof, value);
                }
            }
        }
    }
}



//! @brief adds \f$(\boldsymbol{M}_{JJ} - \boldsymbol{C}_{mat}^T\,\boldsymbol{M}_{KJ} - \boldsymbol{M}_{JK}\,\boldsymbol{C}_{mat} + \boldsymbol{C}_{mat}^T\,\boldsymbol{M}_{KK}\,\boldsymbol{C}_{mat})\,c\f$ to rHessian
//! @remark only calculates active dof types
//! @param rHessian ... global hessian
//! @param rCmat ... constraint matrix
//! @param rScalar ... option to scale the terms
void NuTo::StructureOutputBlockMatrix::ApplyCMatrixScal(BlockSparseMatrix& rHessian, const BlockSparseMatrix& rCmat, double rScalar) const
{

    rHessian.AddScal(JJ, rScalar);

    if (not JJ.GetDofStatus().HasInteractingConstraints())
        return;

    const auto& activeDofTypes = JJ.GetDofStatus().GetActiveDofTypes();

    for (auto i : activeDofTypes)
    {
        for (auto j : activeDofTypes)
        {
            rHessian(i,j).Sub_TransA_B_Plus_C_D_Scal( rCmat(i,i), KJ(i,j), JK(i,j), rCmat(j,j), rScalar);
            rHessian(i,j).Add_TransA_B_C_Scal(        rCmat(i,i), KK(i,j),          rCmat(j,j), rScalar);
        }
    }
}

void NuTo::StructureOutputBlockMatrix::ApplyCMatrix(const BlockSparseMatrix& rCmat)
{

    if (not JJ.GetDofStatus().HasInteractingConstraints())
        return;

    const auto& activeDofTypes = JJ.GetDofStatus().GetActiveDofTypes();

    for (auto i : activeDofTypes)
    {
        for (auto j : activeDofTypes)
        {
            JJ(i,j).Sub_TransA_B_Plus_C_D_Scal( rCmat(i,i), KJ(i,j), JK(i,j), rCmat(j,j), 1);
            JJ(i,j).Add_TransA_B_C_Scal(        rCmat(i,i), KK(i,j),          rCmat(j,j), 1);
        }
    }
}


void NuTo::StructureOutputBlockMatrix::Resize(const std::map<Node::eDof, int>& rNumActiveDofsMap, const std::map<Node::eDof, int>& rNumDependentDofsMap)
{
    assert(rNumActiveDofsMap.size() == rNumDependentDofsMap.size());

    JJ.Resize(rNumActiveDofsMap,    rNumActiveDofsMap   );
    JK.Resize(rNumActiveDofsMap,    rNumDependentDofsMap);
    KJ.Resize(rNumDependentDofsMap, rNumActiveDofsMap   );
    KK.Resize(rNumDependentDofsMap, rNumDependentDofsMap);
}


NuTo::StructureOutputBlockVector NuTo::StructureOutputBlockMatrix::operator *(const StructureOutputBlockVector& rRhs) const
{
    StructureOutputBlockVector result(rRhs.J.GetDofStatus(), true);

    result.J = JJ * rRhs.J + JK * rRhs.K;

    if (JJ.GetDofStatus().HasInteractingConstraints())
    {
        result.K = KJ * rRhs.J + KK * rRhs.K;
    }

    return result;
}

void NuTo::StructureOutputBlockMatrix::CheckDimensions() const
{
    JJ.CheckDimensions();
    JK.CheckDimensions();
    KJ.CheckDimensions();
    KK.CheckDimensions();

    if (!(JJ.GetNumRows()    == JK.GetNumRows()    &&
          KJ.GetNumRows()    == KK.GetNumRows()    &&
          JJ.GetNumColumns() == KJ.GetNumColumns() &&
          JK.GetNumColumns() == KK.GetNumColumns()))
    {
        throw NuTo::MechanicsException(std::string( "[") + __PRETTY_FUNCTION__ + "] Mismatch in Block dimensions! \n\n" +
                                                    "| JJ | JK | \n" +
                                                    "|----+----| \n" +
                                                    "| KJ | KK | \n \n" +
                                                    "Submatrix | Rows / Cols \n" +
                                                    "----------------------- \n" +
                                                    "JJ        | " + std::to_string(JJ.GetNumRows()) + " / " + std::to_string(JJ.GetNumColumns()) + "\n"
                                                    "JK        | " + std::to_string(JK.GetNumRows()) + " / " + std::to_string(JK.GetNumColumns()) + "\n"
                                                    "KJ        | " + std::to_string(KJ.GetNumRows()) + " / " + std::to_string(KJ.GetNumColumns()) + "\n"
                                                    "KK        | " + std::to_string(KK.GetNumRows()) + " / " + std::to_string(KK.GetNumColumns()) + "\n");
    }
}


void NuTo::StructureOutputBlockMatrix::AddScal(const StructureOutputBlockMatrix& rOther, double rFactor)
{
    JJ.AddScal(rOther.JJ, rFactor);
    JK.AddScal(rOther.JK, rFactor);

    if (JJ.GetDofStatus().HasInteractingConstraints())
    {
        KJ.AddScal(rOther.KJ, rFactor);
        KK.AddScal(rOther.KK, rFactor);
    }
}

void NuTo::StructureOutputBlockMatrix::AddScalDiag(const StructureOutputBlockVector& rOther, double rFactor)
{
    JJ.AddScalDiag(rOther.J, rFactor);
    if (JJ.GetDofStatus().HasInteractingConstraints())
    {
        KK.AddScalDiag(rOther.K, rFactor);
    }
}

void NuTo::StructureOutputBlockMatrix::SetZero()
{
    JJ.SetZero();
    JK.SetZero();

    if (JJ.GetDofStatus().HasInteractingConstraints())
    {
        KJ.SetZero();
        KK.SetZero();
    }
}

void NuTo::StructureOutputBlockMatrix::CwiseInvert()
{
    JJ.CwiseInvert();
    JK.CwiseInvert();
    KJ.CwiseInvert();
    KK.CwiseInvert();
}

namespace NuTo
{
std::ostream& operator<<(std::ostream &rOut, const NuTo::StructureOutputBlockMatrix &rStructureOutputBlockMatrix)
{
    rOut << "Active Dofs --- Active Dofs" << std::endl;
    rOut << rStructureOutputBlockMatrix.JJ << std::endl;
    rOut << "Active Dofs --- Dependent Dofs" << std::endl;
    rOut << rStructureOutputBlockMatrix.JK << std::endl;
    rOut << "Dependent Dofs --- Active Dofs" << std::endl;
    rOut << rStructureOutputBlockMatrix.KJ << std::endl;
    rOut << "Dependent Dofs --- Dependent Dofs" << std::endl;
    rOut << rStructureOutputBlockMatrix.KK << std::endl;
    return rOut;
}
}
