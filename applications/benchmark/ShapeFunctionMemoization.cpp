#include "Benchmark.h"

#include "mechanics/integrationtypes/IntegrationType3D4NGauss4Ip.h"
#include "mechanics/interpolationtypes/Interpolation3DTetrahedron.h"
#include "mechanics/nodes/NodeEnum.h"
#include "mechanics/interpolationtypes/InterpolationTypeEnum.h"

#include "math/NaturalCoordinateMemoizer.h"

using Interpolation = NuTo::Interpolation3DTetrahedron;
using Integration = NuTo::IntegrationType3D4NGauss4Ip;
auto dof = NuTo::Node::eDof::DISPLACEMENTS;

BENCHMARK(Hash, WithoutMemoization, runner)
{
    Integration integration;
    Interpolation interpolation(dof, NuTo::Interpolation::eTypeOrder::EQUIDISTANT2, 3);

    while (runner.KeepRunningTime(1))
    {
        for (int i = 0; i < integration.GetNumIntegrationPoints(); ++i)
        {
            auto coords = integration.GetLocalIntegrationPointCoordinates(i);
            auto N = interpolation.CalculateMatrixN(coords);
        }
    }
}

//! @brief current NuTo implementation
BENCHMARK(Hash, IdVectorHash, runner)
{
    Integration integration;
    Interpolation interpolation(dof, NuTo::Interpolation::eTypeOrder::EQUIDISTANT2, 3);
    interpolation.UpdateIntegrationType(integration);

    while (runner.KeepRunningTime(1))
    {
        for (int i = 0; i < integration.GetNumIntegrationPoints(); ++i)
        {
            const auto& N = interpolation.GetMatrixN(i);
            (void)N;
        }
    }
}


BENCHMARK(Hash, NaturalCoordinateHash, runner)
{
    Integration integration;
    Interpolation interpolation(dof, NuTo::Interpolation::eTypeOrder::EQUIDISTANT2, 3);

    std::function<Eigen::MatrixXd(Eigen::Vector3d)> fct = [=](Eigen::Vector3d v) {
        return interpolation.CalculateMatrixN(v);
    };
    NuTo::NaturalCoordinateMemoizer<Eigen::MatrixXd, Eigen::VectorXd> myhash(fct);

    while (runner.KeepRunningTime(1))
    {
        for (int i = 0; i < integration.GetNumIntegrationPoints(); ++i)
        {
            auto coords = integration.GetLocalIntegrationPointCoordinates(i);
            const auto& N = myhash.Get(coords);
            (void)N;
        }
    }
}
