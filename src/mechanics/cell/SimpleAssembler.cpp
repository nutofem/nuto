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

GlobalDofVector ProperlyResizedGlobalVector(DofInfo dofInfo, std::vector<DofType> dofTypes)
{
    GlobalDofVector v;
    for (auto dof : dofTypes)
    {
        v.J[dof].setZero(dofInfo.numIndependentDofs[dof]);
        v.K[dof].setZero(dofInfo.numDependentDofs[dof]);
    }
    return v;
}

GlobalDofMatrixSparse ProperlyResizedGlobalMatrix(DofInfo dofInfo, std::vector<DofType> dofTypes)
{
    GlobalDofMatrixSparse m;
    for (auto dofI : dofTypes)
        for (auto dofJ : dofTypes)
        {
            m.JJ(dofI, dofJ).resize(dofInfo.numIndependentDofs[dofI], dofInfo.numIndependentDofs[dofJ]);
            m.JK(dofI, dofJ).resize(dofInfo.numIndependentDofs[dofI], dofInfo.numDependentDofs[dofJ]);
            m.KJ(dofI, dofJ).resize(dofInfo.numDependentDofs[dofI], dofInfo.numIndependentDofs[dofJ]);
            m.KK(dofI, dofJ).resize(dofInfo.numDependentDofs[dofI], dofInfo.numDependentDofs[dofJ]);
        }
    return m;
}

void ThrowOnZeroDofNumbering(DofInfo dofInfo, std::vector<DofType> dofTypes)
{
    for (auto dof : dofTypes)
        if (not dofInfo.numIndependentDofs.Has(dof))
            throw Exception("[NuTo::SimpleAssembler]",
                            "You did not provide a dof numbering for DofType " + dof.GetName() +
                                    ". Please do so by SimpleAssembler::calling SetDofInfo(...).");
}

GlobalDofVector SimpleAssembler::BuildVector(DofVectorGenerator& entries) const
{
    auto dofTypes = entries.Dofs();
    ThrowOnZeroDofNumbering(mDofInfo, dofTypes);

    GlobalDofVector gradient = ProperlyResizedGlobalVector(mDofInfo, dofTypes);
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
    ThrowOnZeroDofNumbering(mDofInfo, dofTypes);

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
    GlobalDofMatrixSparse hessian = ProperlyResizedGlobalMatrix(mDofInfo, dofTypes);
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
