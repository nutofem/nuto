//
// Created by Thomas Titscher on 03/01/17.
//

#include "Benchmark.h"
#include <eigen3/Eigen/Core>
#include <eigen3/Eigen/Sparse>
#include <random>
#include "math/SparseMatrixCSRVector2General.h"
#include "math/SparseMatrixCSRGeneral.h"

constexpr int numElements = 10000;
constexpr int numDofsPerElement = 3 * 20; // Brick order 2

constexpr int numDofs = 139623; // approx from LinearElasticity test

std::uniform_int_distribution<> distribution(0, numDofs-1);
std::mt19937 generator;


//! @brief dummy element matrix
Eigen::MatrixXd Ke = Eigen::MatrixXd::Random(numDofsPerElement,numDofsPerElement);

//! @brief dummy element dof numbering
Eigen::VectorXi GetGlobalDofNumbers()
{
    Eigen::VectorXi random(numDofsPerElement);
    for (int i = 0; i < random.rows(); ++i)
        random[i] = distribution(generator);
    return random;
}

//! @brief get an idea for the cost of GetGlobalDofNumbers
//! @remark ... its not important
BENCHMARK(AssembleSparseMatrix, GetDofNumbers, runner)
{
    while (runner.KeepRunningIterations(1))
    {
        for (int iElement = 0; iElement < numElements; ++iElement)
            GetGlobalDofNumbers();
    }
}


BENCHMARK(AssembleSparseMatrix, NuToVector2, runner)
{
    NuTo::SparseMatrixCSRVector2General<double> mVector2(numDofs, numDofs);
    while (runner.KeepRunningIterations(1))
    {
        for (int iElement = 0; iElement < numElements; ++iElement)
        {
            Eigen::VectorXi dofs = GetGlobalDofNumbers();
            for (int i = 0; i < Ke.rows(); ++i)
                for (int j = 0; j < Ke.rows(); ++j)
                    mVector2.AddValue(dofs[i], dofs[j], Ke(i,j));
        }
        // convert to CSR format
        NuTo::SparseMatrixCSRGeneral<double> m(mVector2);
    }

    int numDoublesNuTo = mVector2.GetNumEntries();
    int numIntsNuTo    = mVector2.GetNumEntries();
    int numVectorsNuTo = 2 * mVector2.GetNumRows() + 2;

    int sizeNuTo
        = sizeof(double) * numDoublesNuTo
            + sizeof(int)    * numIntsNuTo
            + sizeof(std::vector<char*>) * numVectorsNuTo;

    std::cout << " ##############   Approx size NuToVector2:        " << sizeNuTo / 1024 / 1024 << " MB" << std::endl;
}

BENCHMARK(AssembleSparseMatrix, EigenSparse, runner)
{
    std::vector<Eigen::Triplet<double>> coeffs;
    while (runner.KeepRunningIterations(1))
    {
        size_t numEntriesApprox = numElements*numDofsPerElement*numDofsPerElement;
        coeffs.reserve(numEntriesApprox);

        for (int iElement = 0; iElement < numElements; ++iElement)
        {
            Eigen::VectorXi dofs = GetGlobalDofNumbers();
            for (int i = 0; i < Ke.rows(); ++i)
                for (int j = 0; j < Ke.rows(); ++j)
                    coeffs.push_back(Eigen::Triplet<double>(dofs[i], dofs[j], Ke(i,j)));
        }
        Eigen::SparseMatrix<double> m(numDofs, numDofs);
        m.setFromTriplets(coeffs.begin(), coeffs.end());
        m.makeCompressed();
    }
    int sizeofTriplet = sizeof(Eigen::Triplet<double>);
    int sizeEigen = sizeof(std::vector<char*>) + coeffs.size() * sizeofTriplet;

    std::cout << " ##############   Approx size EigenTripletList:   " << sizeEigen / 1024 / 1024 << " MB"  << std::endl;
}
