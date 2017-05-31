#include "Benchmark.h"

#include "mechanics/integrationtypes/IntegrationType3D4NGauss4Ip.h"
#include "mechanics/interpolationtypes/Interpolation3DTetrahedron.h"
#include "mechanics/nodes/NodeEnum.h"
#include "mechanics/interpolationtypes/InterpolationTypeEnum.h"

using Interpolation = NuTo::Interpolation3DTetrahedron;
using Integration = NuTo::IntegrationType3D4NGauss4Ip;
auto dof = NuTo::Node::eDof::DISPLACEMENTS;

//! @brief current NuTo implementation
BENCHMARK(Hash, IdVectorHash, runner)
{
    Integration integration;
    Interpolation interpolation(dof, NuTo::Interpolation::eTypeOrder::EQUIDISTANT2, 3);

    while (runner.KeepRunningTime(1))
    {
        for (int i = 0; i < integration.GetNumIntegrationPoints(); ++i)
        {
            auto coords = integration.GetLocalIntegrationPointCoordinates(i);
            const auto& N = interpolation.MatrixN(coords);
            (void)N;
        }
    }
}

