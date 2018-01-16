#include <benchmark/benchmark.h>
#include <vector>
#include <numeric>
#include "base/Group.h"
#include "base/ContainerView.h"

/*
 */

std::vector<int> V(int n)
{
    std::vector<int> v(n);
    std::iota(v.begin(), v.end(), 0);
    return v;
}

template <typename TContainer>
void LightWork(const TContainer& v)
{
    benchmark::DoNotOptimize(std::accumulate(v.begin(), v.end(), 0));
}

static void VectorLight(benchmark::State& state)
{
    std::vector<int> v = V(state.range(0));
    for (auto _ : state)
    {
        LightWork(v);
    }
    state.SetComplexityN(state.range(0));
}
static void ViewLight(benchmark::State& state)
{
    std::vector<int> v = V(state.range(0));
    for (auto _ : state)
    {
        NuTo::ContainerView<int> view(v);
        LightWork(view);
    }
    state.SetComplexityN(state.range(0));
}
static void GroupLight(benchmark::State& state)
{
    std::vector<int> v = V(state.range(0));
    for (auto _ : state)
    {
        NuTo::Group<int> group;
        for (auto& val : v)
            group.Add(val);
        LightWork(group);
    }
    state.SetComplexityN(state.range(0));
}
BENCHMARK(VectorLight)->Range(2 << 8, 2 << 16)->Complexity();
BENCHMARK(ViewLight)->Range(2 << 8, 2 << 16)->Complexity();
BENCHMARK(GroupLight)->Range(2 << 8, 2 << 16)->Complexity();

template <typename TContainer>
void HeavyWork(const TContainer& v)
{
    double result = 0;
    for (int val : v)
        benchmark::DoNotOptimize(result += std::exp(std::exp(std::sin(val))));
}

static void VectorHeavy(benchmark::State& state)
{
    std::vector<int> v = V(state.range(0));
    for (auto _ : state)
    {
        HeavyWork(v);
    }
    state.SetComplexityN(state.range(0));
}
static void ViewHeavy(benchmark::State& state)
{
    std::vector<int> v = V(state.range(0));
    for (auto _ : state)
    {
        NuTo::ContainerView<int> view(v);
        HeavyWork(view);
    }
    state.SetComplexityN(state.range(0));
}
static void GroupHeavy(benchmark::State& state)
{
    std::vector<int> v = V(state.range(0));
    for (auto _ : state)
    {
        NuTo::Group<int> group;
        for (auto& val : v)
            group.Add(val);
        HeavyWork(group);
    }
    state.SetComplexityN(state.range(0));
}

BENCHMARK(VectorHeavy)->Range(2 << 8, 2 << 16)->Complexity();
BENCHMARK(ViewHeavy)->Range(2 << 8, 2 << 16)->Complexity();
BENCHMARK(GroupHeavy)->Range(2 << 8, 2 << 16)->Complexity();
BENCHMARK_MAIN();
