#pragma once

#include "mechanics/cell/CellInterface.h"
#include "mechanics/nodes/DofMatrixSparse.h"

namespace NuTo
{
class SimpleAssembler
{
public:
    SimpleAssembler(const NuTo::DofContainer<int>& numDofs,
                    const NuTo::DofContainer<std::vector<int>>& constrainedDofNumbers)
        : mNumDofs(numDofs)
        , mConstrainedToZeroDofNumbers(constrainedDofNumbers)
    {
    }

    DofVector<double> BuildGradient(const std::vector<NuTo::CellInterface*>& cells,
                                    const std::vector<NuTo::DofType*>& dofTypes) const
    {
        DofVector<double> gradient = ProperlyResizedGlobalVector();
        for (NuTo::CellInterface* cell : cells)
        {
            const DofVector<int> numbering = cell->DofNumbering();
            const DofVector<double> cellGradient = cell->Gradient();

            for (const DofType* dof : dofTypes)
            {
                const Eigen::VectorXi& numberingDof = numbering[*dof];
                const Eigen::VectorXd& cellGradientDof = cellGradient[*dof];
                Eigen::VectorXd& globalGradientDof = gradient[*dof];
                for (int i = 0; i < numberingDof.rows(); ++i)
                {
                    int globalDofNumber = numberingDof[i];
                    double globalDofValue = cellGradientDof[i];
                    globalGradientDof[globalDofNumber] = globalDofValue;
                }
            }
        }
        return gradient;
    }

    DofMatrixSparse<double> BuildHessian0(const std::vector<NuTo::CellInterface*>& cells,
                                          const std::vector<NuTo::DofType*>& dofTypes) const
    {
        DofMatrixSparse<double> hessian = ProperlyResizedGlobalMatrix();

        for (NuTo::CellInterface* cell : cells)
        {
            const DofVector<int> numbering = cell->DofNumbering();
            const DofMatrix<double> cellHessian = cell->Hessian0();

            for (const DofType* dofI : dofTypes)
            {
                for (const DofType* dofJ : dofTypes)
                {
                    const Eigen::VectorXi& numberingDofI = numbering[*dofI];
                    const Eigen::VectorXi& numberingDofJ = numbering[*dofJ];
                    const Eigen::MatrixXd& cellHessianDof = cellHessian(*dofI, *dofJ);
                    Eigen::SparseMatrix<double>& globalHessianDof = hessian(*dofI, *dofJ);
                    for (int i = 0; i < numberingDofI.rows(); ++i)
                    {
                        for (int j = 0; j < numberingDofJ.rows(); ++j)
                        {
                            int globalDofNumberI = numberingDofI[i];
                            int globalDofNumberJ = numberingDofJ[j];
                            double globalDofValue = cellHessianDof(i, j);
                            globalHessianDof.coeffRef(globalDofNumberI, globalDofNumberJ) = globalDofValue;
                        }
                    }
                }
            }
        }
        return hessian;
    }

private:
    DofVector<double> ProperlyResizedGlobalVector(const std::vector<NuTo::DofType*>& dofTypes) const
    {
        DofVector<double> v;
        for (auto* dof : dofTypes)
            v[*dof].setZero(mNumDofs[*dof]);
        return v;
    }

    DofMatrixSparse<double> ProperlyResizedGlobalMatrix(const std::vector<NuTo::DofType*>& dofTypes) const
    {
        DofMatrixSparse<double> m;
        for (auto* dofI : dofTypes)
            for (auto* dofJ : dofTypes)
                m(*dofI, *dofJ).resize(mNumDofs[*dofI], mNumDofs[*dofJ]);
        return m;
    }

    NuTo::DofContainer<int> mNumDofs;
    NuTo::DofContainer<std::vector<int>> mConstrainedToZeroDofNumbers;
};
} /* NuTo */
