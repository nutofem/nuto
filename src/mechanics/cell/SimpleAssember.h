#pragma once

#include "base/Group.h"
#include "mechanics/cell/CellInterface.h"
#include "mechanics/dofs/DofContainer.h"
#include "mechanics/dofs/GlobalDofVector.h"
#include "mechanics/dofs/GlobalDofMatrixSparse.h"

namespace NuTo
{
class SimpleAssembler
{
public:
    SimpleAssembler(const NuTo::DofContainer<int>& numIndependentDofs, const NuTo::DofContainer<int>& numDependentDofs)
        : mNumIndependentDofs(numIndependentDofs)
        , mNumDependentDofs(numDependentDofs)
    {
    }

    GlobalDofVector BuildVector(const Groups::Group<NuTo::CellInterface>& cells,
                                const std::vector<NuTo::DofType*>& dofTypes, const VectorOperation& op) const
    {
        GlobalDofVector gradient = ProperlyResizedGlobalVector(dofTypes);
        for (NuTo::CellInterface& cell : cells)
        {
            const DofVector<double> cellGradient = cell.Integrate(op);

            for (const DofType* dof : dofTypes)
            {
                Eigen::VectorXi numberingDof = cell.DofNumbering(*dof);
                const Eigen::VectorXd& cellGradientDof = cellGradient[*dof];
                const int numIndependentDofs = mNumIndependentDofs[*dof];
                for (int i = 0; i < numberingDof.rows(); ++i)
                {
                    int globalDofNumber = numberingDof[i];
                    double globalDofValue = cellGradientDof[i];
                    if (globalDofNumber < numIndependentDofs)
                        gradient.J[*dof][globalDofNumber] += globalDofValue;
                    else
                        gradient.K[*dof][globalDofNumber - numIndependentDofs] += globalDofValue;
                }
            }
        }
        return gradient;
    }

    GlobalDofMatrixSparse BuildMatrix(const Groups::Group<NuTo::CellInterface>& cells,
                                      const std::vector<NuTo::DofType*>& dofTypes, const MatrixOperation& op) const
    {
        GlobalDofMatrixSparse hessian = ProperlyResizedGlobalMatrix(dofTypes);

        for (NuTo::CellInterface& cell : cells)
        {
            const DofMatrix<double> cellHessian = cell.Integrate(op);

            for (const DofType* dofI : dofTypes)
            {
                Eigen::VectorXi numberingDofI = cell.DofNumbering(*dofI);
                for (const DofType* dofJ : dofTypes)
                {
                    Eigen::VectorXi numberingDofJ = cell.DofNumbering(*dofJ);
                    const Eigen::MatrixXd& cellHessianDof = cellHessian(*dofI, *dofJ);

                    const int numIndependentDofsI = mNumIndependentDofs[*dofI];
                    const int numIndependentDofsJ = mNumIndependentDofs[*dofJ];

                    for (int i = 0; i < numberingDofI.rows(); ++i)
                    {
                        for (int j = 0; j < numberingDofJ.rows(); ++j)
                        {
                            const int globalDofNumberI = numberingDofI[i];
                            const int globalDofNumberJ = numberingDofJ[j];
                            const double globalDofValue = cellHessianDof(i, j);

                            const bool activeI = globalDofNumberI < numIndependentDofsI;
                            const bool activeJ = globalDofNumberJ < numIndependentDofsJ;

                            if (activeI)
                            {
                                if (activeJ)
                                {
                                    hessian.JJ(*dofI, *dofJ).coeffRef(globalDofNumberI, globalDofNumberJ) +=
                                            globalDofValue;
                                }
                                else
                                {
                                    hessian.JK(*dofI, *dofJ)
                                            .coeffRef(globalDofNumberI, globalDofNumberJ - numIndependentDofsJ) +=
                                            globalDofValue;
                                }
                            }
                            else
                            {
                                if (activeJ)
                                {
                                    hessian.KJ(*dofI, *dofJ)
                                            .coeffRef(globalDofNumberI - numIndependentDofsI, globalDofNumberJ) +=
                                            globalDofValue;
                                }
                                else
                                {
                                    hessian.KK(*dofI, *dofJ)
                                            .coeffRef(globalDofNumberI - numIndependentDofsI,
                                                      globalDofNumberJ - numIndependentDofsJ) += globalDofValue;
                                }
                            } // argh. any better ideas?
                        }
                    }
                }
            }
        }
        return hessian;
    }

private:
    GlobalDofVector ProperlyResizedGlobalVector(const std::vector<NuTo::DofType*>& dofTypes) const
    {
        GlobalDofVector v;
        for (auto* dof : dofTypes)
        {
            v.J[*dof].setZero(mNumIndependentDofs[*dof]);
            v.K[*dof].setZero(mNumDependentDofs[*dof]);
        }
        return v;
    }

    GlobalDofMatrixSparse ProperlyResizedGlobalMatrix(const std::vector<NuTo::DofType*>& dofTypes) const
    {
        GlobalDofMatrixSparse m;
        for (auto* dofI : dofTypes)
            for (auto* dofJ : dofTypes)
            {
                m.JJ(*dofI, *dofJ).resize(mNumIndependentDofs[*dofI], mNumIndependentDofs[*dofJ]);
                m.JK(*dofI, *dofJ).resize(mNumIndependentDofs[*dofI], mNumDependentDofs[*dofJ]);
                m.KJ(*dofI, *dofJ).resize(mNumDependentDofs[*dofI], mNumIndependentDofs[*dofJ]);
                m.KK(*dofI, *dofJ).resize(mNumDependentDofs[*dofI], mNumDependentDofs[*dofJ]);
            }
        return m;
    }

    NuTo::DofContainer<int> mNumIndependentDofs;
    NuTo::DofContainer<int> mNumDependentDofs;
};
} /* NuTo */
