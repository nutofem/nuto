#include "Benchmark.h"

#include "mechanics/dofs/DofVector.h"
#include "mechanics/dofs/DofContainer.h"

class MapVector : public NuTo::DofContainer<Eigen::VectorXd>
{
public:
    MapVector& operator+=(MapVector&& rhs)
    {
        for (const auto& it : rhs.mData)
        {
            auto& thisData = this->mData[it.first];
            if (thisData.rows() == 0)
                thisData = std::move(it.second);
            else
                thisData += it.second;
        }
        return *this;
    }

    friend MapVector operator*(MapVector lhs, double scalar)
    {
        for (auto& it : lhs.mData)
            it.second *= scalar;
        return lhs;
    }
};

NuTo::DofType d0("zero", 2);
NuTo::DofType d1("one", 2);
NuTo::DofType d2("two", 2);

constexpr int numDofs = 60; // e.g. for quadratic brick
constexpr int numIp = 8;
constexpr int numRuns = 10000;

Eigen::VectorXd randomGradient = Eigen::VectorXd::Random(numDofs);


BENCHMARK(BlockTypes, NuToDofVector, runner)
{
    while (runner.KeepRunningIterations(numRuns))
    {
        NuTo::DofVector<double> v;
        for (int i = 0; i < numIp; ++i)
        {
            NuTo::DofVector<double> gradient(3);
            gradient[d0] = randomGradient;
            gradient[d1] = randomGradient;
            gradient[d2] = randomGradient;
            v += std::move(gradient) * 0.3;
        }
    }
}

BENCHMARK(BlockTypes, NuToDofVectorAddScaled, runner)
{
    while (runner.KeepRunningIterations(numRuns))
    {
        NuTo::DofVector<double> v;
        for (int i = 0; i < numIp; ++i)
        {
            NuTo::DofVector<double> gradient(3);
            gradient[d0] = randomGradient;
            gradient[d1] = randomGradient;
            gradient[d2] = randomGradient;
            v.AddScaled(gradient, 0.3);
        }
    }
}
BENCHMARK(BlockTypes, MapVector, runner)
{
    while (runner.KeepRunningIterations(numRuns))
    {
        MapVector v;
        for (int i = 0; i < numIp; ++i)
        {
            MapVector gradient;
            gradient[d0] = randomGradient;
            gradient[d1] = randomGradient;
            gradient[d2] = randomGradient;
            v += gradient * 0.3;
        }
    }
}


BENCHMARK(BlockTypes, EigenVectorPotential, runner)
{
    while (runner.KeepRunningIterations(numRuns))
    {
        Eigen::VectorXd v;
        for (int i = 0; i < numIp; ++i)
        {
            Eigen::VectorXd gradient;
            gradient.conservativeResize(numDofs);
            gradient.segment(0, numDofs) = randomGradient;
            gradient.conservativeResize(2 * numDofs);
            gradient.segment(numDofs, numDofs) = randomGradient;
            gradient.conservativeResize(3 * numDofs);
            gradient.segment(2 * numDofs, numDofs) = randomGradient;
            v += gradient * 0.3;
        }
    }
}
