//
// Created by Thomas Titscher on 03/01/17.
//

#include "Benchmark.h"
#include <Eigen/Core>
#include <Eigen/Sparse>
#include "math/SparseMatrixCSRVector2General.h"
#include "math/SparseMatrixCSRGeneral.h"

#include "LinearElasticBenchmarkStructure.h"

using namespace NuTo;

constexpr int numBenchmarkRuns = 10;

class Setup
{
public:
    Setup()
        : mBenchmarkStructure({10, 10, 100})
        , mS(mBenchmarkStructure.GetStructure())
    {
        mS.NodeBuildGlobalDofs();

        mDofNumbers.resize(NumElements());
        for (int i = 0; i < NumElements(); ++i)
            mDofNumbers[i] = mS.ElementBuildGlobalDofsColumn(i)[Node::eDof::DISPLACEMENTS];

        mNumDofsPerElement = mDofNumbers[0].rows();
        mKe = Eigen::MatrixXd::Random(mNumDofsPerElement, mNumDofsPerElement);
    }

    int NumDofs() const
    {
        return mS.GetNumTotalDofs();
    }
    int NumElements() const
    {
        return mS.GetNumElements();
    }
    Eigen::VectorXi DofNumbering(int rElementId) const
    {
        return mDofNumbers[rElementId];
    }
    Eigen::MatrixXd Ke()
    {
        return mKe;
    }

private:
    NuTo::Benchmark::LinearElasticBenchmarkStructure mBenchmarkStructure;
    Structure& mS;

    std::vector<Eigen::VectorXi> mDofNumbers;
    int mNumDofsPerElement;
    Eigen::MatrixXd mKe;
};


BENCHMARK(Assembly, NuToVector2, runner)
{
    Setup s;
    NuTo::SparseMatrixCSRVector2General<double> mVector2;
    while (runner.KeepRunningIterations(numBenchmarkRuns))
    {
        mVector2.Resize(s.NumDofs(), s.NumDofs());
        for (int iElement = 0; iElement < s.NumElements(); ++iElement)
        {
            auto dofs = s.DofNumbering(iElement);
            auto Ke = s.Ke();

            for (int i = 0; i < Ke.rows(); ++i)
                for (int j = 0; j < Ke.rows(); ++j)
                    mVector2.AddValue(dofs[i], dofs[j], Ke(i, j));
        }
        // convert to CSR format
        NuTo::SparseMatrixCSRGeneral<double> m(mVector2);
    }

    int numDoublesNuTo = mVector2.GetNumEntries();
    int numIntsNuTo = mVector2.GetNumEntries();
    int numVectorsNuTo = 2 * mVector2.GetNumRows() + 2;

    int sizeNuTo =
            sizeof(double) * numDoublesNuTo + sizeof(int) * numIntsNuTo + sizeof(std::vector<char*>) * numVectorsNuTo;

    std::cout << " ##############   Approx size NuToVector2:        " << sizeNuTo / 1024 / 1024 << " MB" << std::endl;
    std::cout << " ############## NumEntries/Size "
              << (double)mVector2.GetNumEntries() / mVector2.GetNumRows() / mVector2.GetNumColumns() << std::endl;
}


BENCHMARK(Assembly, EigenSparseTriplet, runner)
{
    Setup s;
    std::vector<Eigen::Triplet<double>> coeffs;
    while (runner.KeepRunningIterations(numBenchmarkRuns))
    {
        coeffs.clear();
        size_t numEntriesApprox = s.NumElements() * s.Ke().size();
        coeffs.reserve(numEntriesApprox);

        for (int iElement = 0; iElement < s.NumElements(); ++iElement)
        {
            auto dofs = s.DofNumbering(iElement);
            auto Ke = s.Ke();

            for (int i = 0; i < Ke.rows(); ++i)
                for (int j = 0; j < Ke.rows(); ++j)
                    coeffs.push_back(Eigen::Triplet<double>(dofs[i], dofs[j], Ke(i, j)));
        }
        Eigen::SparseMatrix<double> m(s.NumDofs(), s.NumDofs());
        m.setFromTriplets(coeffs.begin(), coeffs.end());
        m.makeCompressed();
    }
    int sizeofTriplet = sizeof(Eigen::Triplet<double>);
    int sizeEigen = sizeof(std::vector<char*>) + coeffs.size() * sizeofTriplet;

    std::cout << " ##############   Approx size EigenTripletList:   " << sizeEigen / 1024 / 1024 << " MB" << std::endl;
}

BENCHMARK(Assembly, EigenSparseInsert, runner)
{
    Setup s;
    while (runner.KeepRunningIterations(numBenchmarkRuns))
    {
        Eigen::SparseMatrix<double> m(s.NumDofs(), s.NumDofs());
        m.reserve(Eigen::VectorXi::Constant(s.NumDofs(), s.Ke().rows() * 20));
        {
            for (int iElement = 0; iElement < s.NumElements(); ++iElement)
            {
                auto dofs = s.DofNumbering(iElement);
                auto Ke = s.Ke();

                for (int i = 0; i < Ke.rows(); ++i)
                    for (int j = 0; j < Ke.rows(); ++j)
                        m.insert(dofs[i], dofs[j]) = Ke(i, j);
            }
        }
        m.makeCompressed();
    }
}
