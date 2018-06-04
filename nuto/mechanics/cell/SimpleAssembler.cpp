#include "nuto/mechanics/cell/SimpleAssembler.h"
#include "nuto/base/Exception.h"

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

std::vector<DofType> DofIntersection(std::vector<DofType> one, std::vector<DofType> two)
{
    std::vector<DofType> intersection;
    for (DofType d1 : one)
        for (DofType d2 : two)
            if (d1.Id() == d2.Id())
            {
                intersection.push_back(d1);
                continue;
            }
    return intersection;
}

DofVector<double> SimpleAssembler::BuildVector(const Group<CellInterface>& cells, std::vector<DofType> dofTypes,
                                               CellInterface::VectorFunction f) const
{
    ThrowOnZeroDofNumbering(dofTypes);

    DofVector<double> gradient = ProperlyResizedVector(dofTypes);

#pragma omp parallel
    {
        DofVector<double> threadlocalgradient = ProperlyResizedVector(dofTypes);
#pragma omp for nowait
        for (auto cellit = cells.begin(); cellit < cells.end(); cellit++)
        {
            const DofVector<double> cellGradient = cellit->Integrate(f);

            for (DofType dof : DofIntersection(cellGradient.DofTypes(), dofTypes))
            {
                Eigen::VectorXi numberingDof = cellit->DofNumbering(dof);
                const Eigen::VectorXd& cellGradientDof = cellGradient[dof];
                for (int i = 0; i < numberingDof.rows(); ++i)
                    threadlocalgradient[dof](numberingDof[i]) += cellGradientDof[i];
            }
        }
#pragma omp critical
        gradient += threadlocalgradient;
    }
    return gradient;
}

DofVector<double> SimpleAssembler::BuildDiagonallyLumpedMatrix(const Group<CellInterface>& cells,
                                                               std::vector<DofType> dofTypes,
                                                               CellInterface::MatrixFunction f) const
{
    DofVector<double> lumpedMatrix = ProperlyResizedVector(dofTypes);

#pragma omp parallel
    {
        DofVector<double> threadlocallumpedMatrix = ProperlyResizedVector(dofTypes);
#pragma omp for nowait
        for (auto cellit = cells.begin(); cellit < cells.end(); cellit++)
        {
            const DofMatrix<double> localMatrix = cellit->Integrate(f);

            for (DofType dof : DofIntersection(localMatrix.DofTypes(), dofTypes))
            {
                Eigen::VectorXd localDiagonalDof = localMatrix(dof, dof).diagonal();
                double diagonalSum = localDiagonalDof.sum();
                double fullSum = localMatrix(dof, dof).sum();
                localDiagonalDof *= (fullSum / diagonalSum) / dof.GetNum();

                Eigen::VectorXi numberingDof = cellit->DofNumbering(dof);

                for (int i = 0; i < numberingDof.rows(); ++i)
                    threadlocallumpedMatrix[dof](numberingDof[i]) += localDiagonalDof[i];
            }
        }
#pragma omp critical
        lumpedMatrix += threadlocallumpedMatrix;
    }
    return lumpedMatrix;
}

DofMatrixSparse<double> SimpleAssembler::BuildMatrix(const Group<CellInterface>& cells, std::vector<DofType> dofTypes,
                                                     CellInterface::MatrixFunction f) const
{
    ThrowOnZeroDofNumbering(dofTypes);

    using TripletList = std::list<Eigen::Triplet<double>>;
    DofMatrixContainer<TripletList> triplets;

#pragma omp parallel
    {
        DofMatrixContainer<TripletList> localtriplets;

#pragma omp for nowait
        for (auto cellit = cells.begin(); cellit < cells.end(); cellit++)
        {
            const DofMatrix<double> cellHessian = cellit->Integrate(f);
            auto dofTypesToAssemble = DofIntersection(cellHessian.DofTypes(), dofTypes);

            for (DofType dofI : dofTypesToAssemble)
            {
                Eigen::VectorXi numberingDofI = cellit->DofNumbering(dofI);
                for (DofType dofJ : dofTypesToAssemble)
                {
                    Eigen::VectorXi numberingDofJ = cellit->DofNumbering(dofJ);
                    const Eigen::MatrixXd& cellHessianDof = cellHessian(dofI, dofJ);

                    for (int i = 0; i < numberingDofI.rows(); ++i)
                    {
                        for (int j = 0; j < numberingDofJ.rows(); ++j)
                        {
                            const int globalDofNumberI = numberingDofI[i];
                            const int globalDofNumberJ = numberingDofJ[j];
                            const double globalDofValue = cellHessianDof(i, j);

                            localtriplets(dofI, dofJ).push_back({globalDofNumberI, globalDofNumberJ, globalDofValue});
                        }
                    }
                }
            }
        }
#pragma omp critical
        {
            for (DofType dofI : dofTypes)
                for (DofType dofJ : dofTypes)
                {
                    triplets(dofI, dofJ).splice(triplets(dofI, dofJ).end(), localtriplets(dofI, dofJ));
                }
        }
    }
    DofMatrixSparse<double> hessian = ProperlyResizedMatrix(dofTypes);
    for (DofType dofI : dofTypes)
        for (DofType dofJ : dofTypes)
        {
            hessian(dofI, dofJ).setFromTriplets(triplets(dofI, dofJ).begin(), triplets(dofI, dofJ).end());
        }
    return hessian;
}

DofVector<double> SimpleAssembler::ProperlyResizedVector(std::vector<DofType> dofTypes) const
{
    DofVector<double> v;
    for (auto dof : dofTypes)
    {
        v[dof].setZero(mDofInfo.numIndependentDofs[dof] + mDofInfo.numDependentDofs[dof]);
    }
    return v;
}

DofMatrixSparse<double> SimpleAssembler::ProperlyResizedMatrix(std::vector<DofType> dofTypes) const
{
    DofMatrixSparse<double> m;
    for (auto dofI : dofTypes)
        for (auto dofJ : dofTypes)
        {
            m(dofI, dofJ)
                    .resize(mDofInfo.numIndependentDofs[dofI] + mDofInfo.numDependentDofs[dofI],
                            mDofInfo.numIndependentDofs[dofJ] + mDofInfo.numDependentDofs[dofJ]);
        }
    return m;
}
