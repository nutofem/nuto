#pragma once

#include "mechanics/cell/CellInterface.h"
#include "mechanics/dofs/GlobalDofVector.h"
#include "mechanics/dofs/GlobalDofMatrixSparse.h"

namespace NuTo
{
class SimpleAssembler
{
public:
    SimpleAssembler(const NuTo::DofContainer<int>& numActiveDofs, const NuTo::DofContainer<int>& numDependentDofs)
        : mNumActiveDofs(numActiveDofs)
        , mNumDependentDofs(numDependentDofs)
    {
    }

    GlobalDofVector BuildGradient(const std::vector<NuTo::CellInterface*>& cells,
                                  const std::vector<NuTo::DofType*>& dofTypes) const
    {
        GlobalDofVector gradient = ProperlyResizedGlobalVector(dofTypes);
        for (NuTo::CellInterface* cell : cells)
        {
            const DofVector<int> numbering = cell->DofNumbering();
            const DofVector<double> cellGradient = cell->Gradient();

            for (const DofType* dof : dofTypes)
            {
                const Eigen::VectorXi& numberingDof = numbering[*dof];
                const Eigen::VectorXd& cellGradientDof = cellGradient[*dof];
                const int numActiveDofs = mNumActiveDofs[*dof];
                for (int i = 0; i < numberingDof.rows(); ++i)
                {
                    int globalDofNumber = numberingDof[i];
                    double globalDofValue = cellGradientDof[i];
                    if (globalDofNumber < numActiveDofs)
                        gradient.J[*dof][globalDofNumber] += globalDofValue;
                    else
                        gradient.K[*dof][globalDofNumber - numActiveDofs] += globalDofValue;
                }
            }
        }
        return gradient;
    }

    GlobalDofMatrixSparse BuildHessian0(const std::vector<NuTo::CellInterface*>& cells,
                                        const std::vector<NuTo::DofType*>& dofTypes) const
    {
        GlobalDofMatrixSparse hessian = ProperlyResizedGlobalMatrix(dofTypes);

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

                    const int numActiveDofsI = mNumActiveDofs[*dofI];
                    const int numActiveDofsJ = mNumActiveDofs[*dofJ];

                    for (int i = 0; i < numberingDofI.rows(); ++i)
                    {
                        for (int j = 0; j < numberingDofJ.rows(); ++j)
                        {
                            const int globalDofNumberI = numberingDofI[i];
                            const int globalDofNumberJ = numberingDofJ[j];
                            const double globalDofValue = cellHessianDof(i, j);

                            const bool activeI = globalDofNumberI < numActiveDofsI;
                            const bool activeJ = globalDofNumberJ < numActiveDofsJ;

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
                                            .coeffRef(globalDofNumberI, globalDofNumberJ - numActiveDofsJ) +=
                                            globalDofValue;
                                }
                            }
                            else
                            {
                                if (activeJ)
                                {
                                    hessian.KJ(*dofI, *dofJ)
                                            .coeffRef(globalDofNumberI - numActiveDofsI, globalDofNumberJ) +=
                                            globalDofValue;
                                }
                                else
                                {
                                    hessian.KK(*dofI, *dofJ)
                                            .coeffRef(globalDofNumberI - numActiveDofsI,
                                                      globalDofNumberJ - numActiveDofsJ) += globalDofValue;
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
            v.J[*dof].setZero(mNumActiveDofs[*dof]);
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
                m.JJ(*dofI, *dofJ).resize(mNumActiveDofs[*dofI], mNumActiveDofs[*dofJ]);
                m.JK(*dofI, *dofJ).resize(mNumActiveDofs[*dofI], mNumDependentDofs[*dofJ]);
                m.KJ(*dofI, *dofJ).resize(mNumDependentDofs[*dofI], mNumActiveDofs[*dofJ]);
                m.KK(*dofI, *dofJ).resize(mNumDependentDofs[*dofI], mNumDependentDofs[*dofJ]);
            }
        return m;
    }

    NuTo::DofContainer<int> mNumActiveDofs;
    NuTo::DofContainer<int> mNumDependentDofs;
};
} /* NuTo */
