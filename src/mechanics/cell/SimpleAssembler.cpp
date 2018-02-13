#include "mechanics/cell/SimpleAssembler.h"

using namespace NuTo;

SimpleAssembler::SimpleAssembler(DofContainer<int> numIndependentDofs, DofContainer<int> numDependentDofs)
    : mNumIndependentDofs(numIndependentDofs)
    , mNumDependentDofs(numDependentDofs)
{
}

GlobalDofVector SimpleAssembler::BuildVector(const Group<CellInterface>& cells, std::vector<DofType> dofTypes,
                                             CellInterface::VectorFunction f) const
{
    GlobalDofVector gradient = ProperlyResizedGlobalVector(dofTypes);
    for (NuTo::CellInterface& cell : cells)
    {
        const DofVector<double> cellGradient = cell.Integrate(f);

        for (DofType dof : dofTypes)
        {
            Eigen::VectorXi numberingDof = cell.DofNumbering(dof);
            const Eigen::VectorXd& cellGradientDof = cellGradient[dof];
            for (int i = 0; i < numberingDof.rows(); ++i)
                gradient(dof, numberingDof[i]) += cellGradientDof[i];
        }
    }
    return gradient;
}

GlobalDofVector SimpleAssembler::BuildDiagonallyLumpedMatrix(const Group<CellInterface>& cells,
                                                             std::vector<DofType> dofTypes,
                                                             CellInterface::MatrixFunction f) const
{
    GlobalDofVector lumpedMatrix = ProperlyResizedGlobalVector(dofTypes);
    for (NuTo::CellInterface& cell : cells)
    {
        const DofMatrix<double> localMatrix = cell.Integrate(f);

        for (DofType dof : dofTypes)
        {
            Eigen::VectorXd localDiagonalDof = localMatrix(dof, dof).diagonal();
            double diagonalSum = localDiagonalDof.sum();
            double fullSum = localMatrix(dof, dof).sum();
            localDiagonalDof *= (fullSum / diagonalSum) / dof.GetNum();

            Eigen::VectorXi numberingDof = cell.DofNumbering(dof);

            for (int i = 0; i < numberingDof.rows(); ++i)
                lumpedMatrix(dof, numberingDof[i]) += localDiagonalDof[i];
        }
    }
    return lumpedMatrix;
}

GlobalDofMatrixSparse SimpleAssembler::BuildMatrix(const Group<CellInterface>& cells, std::vector<DofType> dofTypes,
                                                   CellInterface::MatrixFunction f) const
{
    using TripletList = std::vector<Eigen::Triplet<double>>;
    DofMatrixContainer<TripletList> tripletsJJ;
    DofMatrixContainer<TripletList> tripletsJK;
    DofMatrixContainer<TripletList> tripletsKJ;
    DofMatrixContainer<TripletList> tripletsKK;

    for (NuTo::CellInterface& cell : cells)
    {
        const DofMatrix<double> cellHessian = cell.Integrate(f);

        for (DofType dofI : dofTypes)
        {
            Eigen::VectorXi numberingDofI = cell.DofNumbering(dofI);
            for (DofType dofJ : dofTypes)
            {
                Eigen::VectorXi numberingDofJ = cell.DofNumbering(dofJ);
                const Eigen::MatrixXd& cellHessianDof = cellHessian(dofI, dofJ);

                const int numIndependentDofsI = mNumIndependentDofs[dofI];
                const int numIndependentDofsJ = mNumIndependentDofs[dofJ];

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
        v.J[dof].setZero(mNumIndependentDofs[dof]);
        v.K[dof].setZero(mNumDependentDofs[dof]);
    }
    return v;
}

GlobalDofMatrixSparse SimpleAssembler::ProperlyResizedGlobalMatrix(std::vector<DofType> dofTypes) const
{
    GlobalDofMatrixSparse m;
    for (auto dofI : dofTypes)
        for (auto dofJ : dofTypes)
        {
            m.JJ(dofI, dofJ).resize(mNumIndependentDofs[dofI], mNumIndependentDofs[dofJ]);
            m.JK(dofI, dofJ).resize(mNumIndependentDofs[dofI], mNumDependentDofs[dofJ]);
            m.KJ(dofI, dofJ).resize(mNumDependentDofs[dofI], mNumIndependentDofs[dofJ]);
            m.KK(dofI, dofJ).resize(mNumDependentDofs[dofI], mNumDependentDofs[dofJ]);
        }
    return m;
}
