#include <benchmark/benchmark.h>
#include "mechanics/constraints/Constraints.h"

/*
 * The constraint creation requires checks with all previously created constraints. In a naive implementation, this has
 * O(N^2) complexity, which is bad (too slow). This benchmark calculates the complexity - hopefully to be O(N log N).
 */

auto zero = [](double) { return 0.; };
const NuTo::DofType d("...", 1);

void Constraints(benchmark::State& state)
{
    std::vector<NuTo::NodeSimple> nodes(state.range(0), NuTo::NodeSimple(0));
    for (int i = 0; i < nodes.size(); ++i)
        nodes[i].SetDofNumber(0, i);

    for (auto _ : state)
    {
        NuTo::Constraint::Constraints c;
        for (int i = 0; i < nodes.size(); ++i)
            c.Add(d, {nodes[i], 0, zero});
    }
    state.SetComplexityN(state.range(0));
}

BENCHMARK(Constraints)->RangeMultiplier(10)->Range(1, 1e6)->Complexity();
BENCHMARK_MAIN()
