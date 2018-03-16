#include <benchmark/benchmark.h>
#include "nuto/mechanics/dofs/DofMatrixSparseConvertEigen.h"
#include <iostream>

using namespace NuTo;

struct TestMatrix
{
    DofType dof0 = DofType("dof0", 1);
    DofType dof1 = DofType("dof1", 1);
    DofMatrixSparse<double> m;

    TestMatrix(int n, int width)
    {
        m(dof0, dof0) = BandMatrix(n, width);
        m(dof0, dof1) = BandMatrix(n, width);
        m(dof1, dof0) = BandMatrix(n, width);
        m(dof1, dof1) = BandMatrix(n, width);
    }

private:
    Eigen::SparseMatrix<double> BandMatrix(int n, int width)
    {
        std::vector<Eigen::Triplet<double>> triplets;
        triplets.reserve(n * ((width - 1) * 2 + 1));
        for (int i = 0; i < n; ++i)
        {
            for (int j = i - width + 1; j < i + width; ++j)
            {
                int row = i;
                int col = std::max(std::min(j, n - 1), 0);
                triplets.emplace_back(Eigen::Triplet<double>(row, col, 42.));
            }
        }
        Eigen::SparseMatrix<double> m(n, n);
        m.setFromTriplets(triplets.begin(), triplets.end());
        m.makeCompressed();
        return m;
    }
};

//! @brief exports two dof types
static void TwoDofs(benchmark::State& s)
{
    TestMatrix t(s.range(0), 5);
    for (auto _ : s)
        auto m = ToEigen(t.m, {t.dof0, t.dof1});

    s.SetComplexityN(s.range(0));
}
BENCHMARK(TwoDofs)->RangeMultiplier(10)->Range(10, 1e6)->Complexity();

//! @brief exports a single dof type.
//! @remark ToEigen avoids expensive export magic, but is forced to make a copy of the matrix. This, prefer the usage in
//! `OneDofFast`
static void OneDof(benchmark::State& s)
{
    TestMatrix t(s.range(0), 5);
    for (auto _ : s)
        benchmark::DoNotOptimize(ToEigen(t.m, {t.dof1}));

    s.SetComplexityN(s.range(0));
}
BENCHMARK(OneDof)->RangeMultiplier(10)->Range(10, 1e6)->Complexity();

//! @brief Just a reminder that there is no need to _export_ a single dof. Access via operator(dof, dof) is blazing
//! fast.
static void OneDofFast(benchmark::State& s)
{
    TestMatrix t(s.range(0), 5);
    for (auto _ : s)
        benchmark::DoNotOptimize(t.m(t.dof1, t.dof1));

    s.SetComplexityN(s.range(0));
}
BENCHMARK(OneDofFast)->RangeMultiplier(10)->Range(10, 1e6)->Complexity();

BENCHMARK_MAIN();
