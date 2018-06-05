#include <iostream>

#include "nuto/mechanics/mesh/MeshCompanion.h"
#include "nuto/mechanics/nodes/NodeSimple.h"
#include "nuto/mechanics/mesh/MeshFem.h"
#include "nuto/mechanics/interpolation/InterpolationTriangleQuadratic.h"

using namespace NuTo;

int main(int, char* argv[])
{
    std::cout << "Create a test mesh" << std::endl;

    MeshFem mesh;

    auto& nd0 = mesh.Nodes.Add(Eigen::Vector2d(0., 0.));
    auto& nd1 = mesh.Nodes.Add(Eigen::Vector2d(1., 0.));
    auto& nd2 = mesh.Nodes.Add(Eigen::Vector2d(1., 1.));
    auto& nd3 = mesh.Nodes.Add(Eigen::Vector2d(0., 1.));

    auto& ndA1 = mesh.Nodes.Add(Eigen::Vector2d(0.5, 0.));
    auto& ndA2 = mesh.Nodes.Add(Eigen::Vector2d(1.0, 0.5));
    auto& ndA3 = mesh.Nodes.Add(Eigen::Vector2d(0.5, 1.0));
    auto& ndA4 = mesh.Nodes.Add(Eigen::Vector2d(0., 0.5));
    auto& ndA5 = mesh.Nodes.Add(Eigen::Vector2d(0.5, 0.5));

    InterpolationTriangleQuadratic ipol;

    mesh.Elements.Add({{{nd0, nd1, nd2, ndA1, ndA2, ndA5}, ipol}});
    mesh.Elements.Add({{{nd2, nd3, nd0, ndA3, ndA4, ndA5}, ipol}});

    Group<ElementCollectionFem> allElements = mesh.ElementsTotal();

    auto edges0 = AddEdgeElements(&mesh, mesh.Elements[0]);
    std::cout << "Added edges of element 0   : " << edges0.Size() << std::endl;

    auto edgesAll = AddEdgeElements(&mesh, allElements);
    std::cout << "Added edges of all elements: " << edgesAll.Size() << std::endl;
}
