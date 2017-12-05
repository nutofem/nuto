#include "Benchmark.h"
#include "mechanics/constraints/Constraints.h"

auto zero = [](double) { return 0.; };
const NuTo::DofType d("...", 1);


void Run(BenchmarkInternal::Runner& runner, int size)
{
    std::vector<NuTo::NodeSimple> nodes(size, NuTo::NodeSimple(0));
    for (int i = 0; i < nodes.size(); ++i)
        nodes[i].SetDofNumber(0, i);

    while (runner.KeepRunningTime(0.1))
    {
        NuTo::Constraint::Constraints c;
        for (int i = 0; i < nodes.size(); ++i)
            c.Add(d, {nodes[i], 0, zero});
    }
}

BENCHMARK(Constraints, c100, runner)
{
    Run(runner, 100);
};

BENCHMARK(Constraints, c1000, runner)
{
    Run(runner, 1000);
};

BENCHMARK(Constraints, c10000, runner)
{
    Run(runner, 10000);
};

BENCHMARK(Constraints, c100000, runner)
{
    Run(runner, 100000);
};

BENCHMARK(Constraints, c1000000, runner)
{
    Run(runner, 1000000);
};
