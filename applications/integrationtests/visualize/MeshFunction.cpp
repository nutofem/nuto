#include <cmath>
#include <boost/ptr_container/ptr_vector.hpp>
#include "nuto/mechanics/mesh/UnitMeshFem.h"
#include "nuto/mechanics/mesh/GeometryMeshFem.h"
#include "nuto/mechanics/mesh/MeshFem.h"
#include "nuto/mechanics/cell/Cell.h"
#include "nuto/mechanics/integrationtypes/IntegrationTypeTensorProduct.h"
#include "nuto/visualize/AverageHandler.h"
#include "nuto/visualize/Visualizer.h"

constexpr double pi()
{
    return std::atan(1) * 4;
}

using namespace NuTo;

int main()
{
    auto meshTmp = UnitMeshFem::CreateQuads(50, 50);
    auto geoMesh = UnitMeshFem::Transform(std::move(meshTmp), [](Eigen::VectorXd coords) { return 2 * pi() * coords; });
    MeshFem mesh(geoMesh);
    boost::ptr_vector<CellInterface> cells;
    IntegrationTypeTensorProduct<2> integrationType(2, eIntegrationMethod::GAUSS);
    int cellId = 0;
    Group<CellInterface> visualizationCells;
    for (size_t i = 0; i < mesh.NumElements(); i++)
    {
        auto& element = mesh.GetElement(i);
        cells.push_back(new Cell(element, integrationType, cellId++));
        visualizationCells.Add(cells.back());
    }

    Visualize::Visualizer visualize(visualizationCells, Visualize::AverageHandler());

    auto sine_function = [](Eigen::VectorXd x) { return Eigen::Vector3d(sin(x[0]), sin(x[1]), sin(x[0] * x[1])); };
    visualize.PointData(sine_function, "Sine");
    visualize.WriteVtuFile("output.vtu");
}
