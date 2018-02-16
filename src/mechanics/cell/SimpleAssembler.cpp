#include "mechanics/cell/SimpleAssembler.h"
#include "base/Exception.h"

using namespace NuTo;

SimpleAssembler::SimpleAssembler(DofInfo dofInfo)
    : mDofInfo(dofInfo)
{
}

void SimpleAssembler::SetDofInfo(DofInfo dofInfo)
{
    mDofInfo = dofInfo;
}

void SimpleAssembler::ThrowOnZeroDofNumbering(std::vector<DofType> dofTypes) const
{
    for (auto dof : dofTypes)
        if (not mDofInfo.numIndependentDofs.Has(dof))
            throw Exception("[NuTo::SimpleAssembler]",
                            "You did not provide a dof numbering for DofType " + dof.GetName() +
                                    ". Please do so by SimpleAssembler::calling SetDofInfo(...).");
}

GlobalDofVector SimpleAssembler::BuildVector(DofVectorGenerator& entries) const
{
    auto dofTypes = entries.Dofs();
    ThrowOnZeroDofNumbering(dofTypes);

    GlobalDofVector gradient = ProperlyResizedGlobalVector(dofTypes);
    for (VectorEntry& entry : entries)
    {
        const DofVector<double>& cellGradient = entry.first;

        for (DofType dof : dofTypes)
        {
            const Eigen::VectorXi& numberingDof = entry.second[dof];
            const Eigen::VectorXd& cellGradientDof = cellGradient[dof];
            for (int i = 0; i < numberingDof.rows(); ++i)
                gradient(dof, numberingDof[i]) += cellGradientDof[i];
        }
    }
    return gradient;
}

GlobalDofMatrixSparse SimpleAssembler::BuildMatrix(DofMatrixGenerator& entries) const
{
    auto dofTypes = entries.Dofs();
    ThrowOnZeroDofNumbering(dofTypes);

    using TripletList = std::vector<Eigen::Triplet<double>>;
    DofMatrixContainer<TripletList> tripletsJJ;
    DofMatrixContainer<TripletList> tripletsJK;
    DofMatrixContainer<TripletList> tripletsKJ;
    DofMatrixContainer<TripletList> tripletsKK;

    for (MatrixEntry& entry : entries)
    {
        const DofMatrix<double>& cellHessian = entry.first;
        ;
        const DofVector<int>& dofNumbering = entry.second;
        ;

        for (DofType dofI : dofTypes)
        {
            const Eigen::VectorXi& numberingDofI = dofNumbering[dofI];
            for (DofType dofJ : dofTypes)
            {
                const Eigen::VectorXi& numberingDofJ = dofNumbering[dofJ];
                const Eigen::MatrixXd& cellHessianDof = cellHessian(dofI, dofJ);

                const int numIndependentDofsI = mDofInfo.numIndependentDofs[dofI];
                const int numIndependentDofsJ = mDofInfo.numIndependentDofs[dofJ];

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
                                tripletsJJ(dofI, dofJ).push_back({globalDofNumberI, globalDofNumberJ, globalDofValue});
                            }
                            else
                            {
                                tripletsJK(dofI, dofJ)
                                        .push_back({globalDofNumberI, globalDofNumberJ - numIndependentDofsJ,
                                                    globalDofValue});
                            }
                        }
                        else
                        {
                            if (activeJ)
                            {
                                tripletsKJ(dofI, dofJ)
                                        .push_back({globalDofNumberI - numIndependentDofsI, globalDofNumberJ,
                                                    globalDofValue});
                            }
                            else
                            {
                                tripletsKK(dofI, dofJ)
                                        .push_back({globalDofNumberI - numIndependentDofsI,
                                                    globalDofNumberJ - numIndependentDofsJ, globalDofValue});
                            }
                        } // argh. any better ideas?
                    }
                }
            }
        }
    }
    GlobalDofMatrixSparse hessian = ProperlyResizedGlobalMatrix(dofTypes);
    for (DofType dofI : dofTypes)
        for (DofType dofJ : dofTypes)
        {
            hessian.JJ(dofI, dofJ).setFromTriplets(tripletsJJ(dofI, dofJ).begin(), tripletsJJ(dofI, dofJ).end());
            hessian.JK(dofI, dofJ).setFromTriplets(tripletsJK(dofI, dofJ).begin(), tripletsJK(dofI, dofJ).end());
            hessian.KJ(dofI, dofJ).setFromTriplets(tripletsKJ(dofI, dofJ).begin(), tripletsKJ(dofI, dofJ).end());
            hessian.KK(dofI, dofJ).setFromTriplets(tripletsKK(dofI, dofJ).begin(), tripletsKK(dofI, dofJ).end());
        }
    return hessian;
}

GlobalDofVector SimpleAssembler::ProperlyResizedGlobalVector(std::vector<DofType> dofTypes) const
{
    GlobalDofVector v;
    for (auto dof : dofTypes)
    {
        v.J[dof].setZero(mDofInfo.numIndependentDofs[dof]);
        v.K[dof].setZero(mDofInfo.numDependentDofs[dof]);
    }
    return v;
}

GlobalDofMatrixSparse SimpleAssembler::ProperlyResizedGlobalMatrix(std::vector<DofType> dofTypes) const
{
    GlobalDofMatrixSparse m;
    for (auto dofI : dofTypes)
        for (auto dofJ : dofTypes)
        {
            m.JJ(dofI, dofJ).resize(mDofInfo.numIndependentDofs[dofI], mDofInfo.numIndependentDofs[dofJ]);
            m.JK(dofI, dofJ).resize(mDofInfo.numIndependentDofs[dofI], mDofInfo.numDependentDofs[dofJ]);
            m.KJ(dofI, dofJ).resize(mDofInfo.numDependentDofs[dofI], mDofInfo.numIndependentDofs[dofJ]);
            m.KK(dofI, dofJ).resize(mDofInfo.numDependentDofs[dofI], mDofInfo.numDependentDofs[dofJ]);
        }
    return m;
}
