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
    for (NuTo::CellInterface& cell : cells)
    {
        const DofVector<double> cellGradient = cell.Integrate(f);

        for (DofType dof : DofIntersection(cellGradient.DofTypes(), dofTypes))
        {
            Eigen::VectorXi numberingDof = cell.DofNumbering(dof);
            const Eigen::VectorXd& cellGradientDof = cellGradient[dof];
            for (int i = 0; i < numberingDof.rows(); ++i)
                gradient[dof](numberingDof[i]) += cellGradientDof[i];
        }
    }
    return gradient;
}

DofVector<double> SimpleAssembler::BuildDiagonallyLumpedMatrix(const Group<CellInterface>& cells,
                                                               std::vector<DofType> dofTypes,
                                                               CellInterface::MatrixFunction f) const
{
    DofVector<double> lumpedMatrix = ProperlyResizedVector(dofTypes);
    for (NuTo::CellInterface& cell : cells)
    {
        const DofMatrix<double> localMatrix = cell.Integrate(f);

        for (DofType dof : DofIntersection(localMatrix.DofTypes(), dofTypes))
        {
            Eigen::VectorXd localDiagonalDof = localMatrix(dof, dof).diagonal();
            double diagonalSum = localDiagonalDof.sum();
            double fullSum = localMatrix(dof, dof).sum();
            localDiagonalDof *= (fullSum / diagonalSum) / dof.GetNum();

            Eigen::VectorXi numberingDof = cell.DofNumbering(dof);

            for (int i = 0; i < numberingDof.rows(); ++i)
                lumpedMatrix[dof](numberingDof[i]) += localDiagonalDof[i];
        }
    }
    return lumpedMatrix;
}

class DofNumberMapping
{
public:
    DofNumberMapping(std::vector<DofType> activeDofTypes, DofInfo dofInfo)
    {
        int offSet = 0;
        mOffSets.resize(activeDofTypes.size());
        auto numAllDofs = dofInfo.numAllDofs();
        for(auto dof : activeDofTypes)
        {
            mOffSets[dof.Id()] = offSet;
            offSet += numAllDofs[dof];
        }
        overallSize = offSet;
    }

    int map(int localId, DofType dof)
    {
        return mOffSets[dof.Id()] + localId; 
    }

    int overallSize;
private:
    std::vector<int> mOffSets;
};

Eigen::SparseMatrix<double> SimpleAssembler::BuildMatrix(const Group<CellInterface>& cells,
                                                         std::vector<DofType> dofTypes,
                                                         DofInfo dofInfo,
                                                         CellInterface::MatrixFunction f) const
{
    ThrowOnZeroDofNumbering(dofTypes);

    using TripletList = std::vector<Eigen::Triplet<double>>;
    TripletList triplets;

    Eigen::SparseMatrix<double> hessian;

    for (NuTo::CellInterface& cell : cells)
    {
        const DofMatrix<double> cellHessian = cell.Integrate(f);
        auto dofTypesToAssemble = DofIntersection(cellHessian.DofTypes(), dofTypes);

        DofNumberMapping mapping(dofTypesToAssemble, dofInfo);
        hessian.resize(mapping.overallSize, mapping.overallSize);
        for (DofType dofI : dofTypesToAssemble)
        {
            Eigen::VectorXi numberingDofI = cell.DofNumbering(dofI);
            for (DofType dofJ : dofTypesToAssemble)
            {
                Eigen::VectorXi numberingDofJ = cell.DofNumbering(dofJ);
                const Eigen::MatrixXd& cellHessianDof = cellHessian(dofI, dofJ);

                for (int i = 0; i < numberingDofI.rows(); ++i)
                {
                    for (int j = 0; j < numberingDofJ.rows(); ++j)
                    {
                        const int globalDofNumberI = mapping.map(numberingDofI[i], dofI);
                        const int globalDofNumberJ = mapping.map(numberingDofJ[j], dofJ);
                        const double globalDofValue = cellHessianDof(i, j);

                        triplets.push_back({globalDofNumberI, globalDofNumberJ, globalDofValue});
                    }
                }
            }
        }
    }

    hessian.setFromTriplets(triplets.begin(), triplets.end());
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
