#include "Benchmark.h"

#include "mechanics/integrationtypes/IntegrationType3D4NGauss4Ip.h"
#include "mechanics/interpolationtypes/Interpolation3DTetrahedron.h"
#include "mechanics/nodes/NodeEnum.h"
#include "mechanics/interpolationtypes/InterpolationTypeEnum.h"
#include "mechanics/elements/ElementShapeFunctions.h"

#include <map>

struct Lobatto
{
    std::vector<Eigen::Vector3d> ips;
    Lobatto()
    {
        std::vector<double> ip1D = {-1., -0.654653670707977087, 0., +0.654653670707977087, +1.};
        int num = ip1D.size();
        ips.reserve(num*num*num);
        for (int i = 0; i < num; i++)
            for (int j = 0; j < num; j++)
                for (int k = 0; k < num; k++)
                    ips.push_back({ip1D[i], ip1D[j], ip1D[k]});
    }
};

auto testFunction = NuTo::ShapeFunctions3D::DerivativeShapeFunctionsBrickSpectralOrder4;
using result = Eigen::Matrix<double, 125, 3>;

template <typename TMemoizer>
void Run(BenchmarkInternal::Runner& runner)
{
    TMemoizer memo(testFunction);
    Lobatto l;
    while (runner.KeepRunningTime(1))
        for (const auto& ip : l.ips)
            memo.Get(ip);
}


BENCHMARK(Hash, Vector, runner)
{
    Run<NuTo::NaturalCoordinateMemoizer<result, Eigen::Vector3d>>(runner);
}

BENCHMARK(Hash, Map, runner)
{
    Run<NuTo::NaturalCoordinateMemoizerMap<result, Eigen::Vector3d>>(runner);
}

BENCHMARK(Hash, Unordered, runner)
{
    Run<NuTo::NaturalCoordinateMemoizerUnorderedMap<result, Eigen::Vector3d>>(runner);
}
