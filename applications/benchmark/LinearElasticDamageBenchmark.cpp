#include <benchmark/benchmark.h>

#include "nuto/mechanics/constitutive/LinearElastic.h"
#include "nuto/mechanics/constitutive/LinearElasticDamage.h"

const double E = 20000;
const double nu = 0.2;
const double omega = 0.2;
const NuTo::EngineeringStrain<3> strain = NuTo::EngineeringStrain<3>::Ones();

static void FullStress(benchmark::State& state)
{
    NuTo::Laws::LinearElastic<3> law(E, nu);
    for (auto _ : state)
        benchmark::DoNotOptimize((1. - omega) * law.Stress(strain));
}
static void FullTangentStrain(benchmark::State& state)
{
    NuTo::Laws::LinearElastic<3> law(E, nu);
    for (auto _ : state)
        benchmark::DoNotOptimize((1. - omega) * law.Tangent(strain));
}
static void FullTangentOmega(benchmark::State& state)
{
    NuTo::Laws::LinearElastic<3> law(E, nu);
    for (auto _ : state)
        benchmark::DoNotOptimize(-law.Tangent(strain));
}


static void UnilateralStress(benchmark::State& state)
{
    NuTo::Laws::LinearElasticDamage<3> law(E, nu);
    for (auto _ : state)
        benchmark::DoNotOptimize(law.Stress(strain, omega));
}

static void UnilateralTangentStrain(benchmark::State& state)
{
    NuTo::Laws::LinearElasticDamage<3> law(E, nu);
    for (auto _ : state)
        benchmark::DoNotOptimize(law.DstressDstrain(strain, omega));
}
static void UnilateralTangentOmega(benchmark::State& state)
{
    NuTo::Laws::LinearElasticDamage<3> law(E, nu);
    for (auto _ : state)
        benchmark::DoNotOptimize(law.DstressDomega(strain, omega));
}


BENCHMARK(FullStress);
BENCHMARK(FullTangentStrain);
BENCHMARK(FullTangentOmega);
BENCHMARK(UnilateralStress);
BENCHMARK(UnilateralTangentStrain);
BENCHMARK(UnilateralTangentOmega);

BENCHMARK_MAIN();
