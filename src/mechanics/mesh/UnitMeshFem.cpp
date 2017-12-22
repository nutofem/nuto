#include "mechanics/mesh/UnitMeshFem.h"
#include "mechanics/interpolation/InterpolationTriangleLinear.h"
#include "mechanics/interpolation/InterpolationTrussLinear.h"
#include "mechanics/interpolation/InterpolationQuadLinear.h"

using namespace NuTo;

// anonymous helper functions for mesh creation and transformation
MeshFem CreateNodes1D(int numX)
{
    MeshFem mesh;
    for (int iX = 0; iX < numX + 1; ++iX)
    {
        const double x = static_cast<double>(iX) / numX;
        mesh.Nodes.Add({x});
    }
    return mesh;
}


MeshFem CreateNodes2D(int numX, int numY)
{
    MeshFem mesh;
    for (int iY = 0; iY < numY + 1; ++iY)
        for (int iX = 0; iX < numX + 1; ++iX)
        {
            const double x = static_cast<double>(iX) / numX;
            const double y = static_cast<double>(iY) / numY;
            mesh.Nodes.Add(Eigen::Vector2d(x, y));
        }
    return mesh;
}


MeshFem UnitMeshFem::CreateLines(int numX)
{
    MeshFem mesh = CreateNodes1D(numX);
    const auto& interpolation = mesh.CreateInterpolation(NuTo::InterpolationTrussLinear());
    for (int i = 0; i < numX; ++i)
    {
        auto& nl = mesh.Nodes[i];
        auto& nr = mesh.Nodes[i + 1];
        mesh.Elements.Add({{{nl, nr}, interpolation}});
    }
    return mesh;
}


MeshFem UnitMeshFem::CreateTriangles(int numX, int numY)
{
    MeshFem mesh = CreateNodes2D(numX, numY);
    const auto& interpolation = mesh.CreateInterpolation(NuTo::InterpolationTriangleLinear());
    for (int iY = 0; iY < numY; ++iY)
        for (int iX = 0; iX < numX; ++iX)
        {
            auto& node0 = mesh.Nodes[iX + iY * (numX + 1)];
            auto& node1 = mesh.Nodes[iX + 1 + iY * (numX + 1)];
            auto& node2 = mesh.Nodes[iX + 1 + (iY + 1) * (numX + 1)];
            auto& node3 = mesh.Nodes[iX + (iY + 1) * (numX + 1)];
            mesh.Elements.Add({{{node0, node1, node2}, interpolation}});
            mesh.Elements.Add({{{node0, node2, node3}, interpolation}});
        }
    return mesh;
}

MeshFem UnitMeshFem::CreateQuads(int numX, int numY)
{
    MeshFem mesh = CreateNodes2D(numX, numY);
    const auto& interpolation = mesh.CreateInterpolation(NuTo::InterpolationQuadLinear());
    for (int iY = 0; iY < numY; ++iY)
        for (int iX = 0; iX < numX; ++iX)
        {
            auto& node0 = mesh.Nodes[iX + iY * (numX + 1)];
            auto& node1 = mesh.Nodes[iX + 1 + iY * (numX + 1)];
            auto& node2 = mesh.Nodes[iX + 1 + (iY + 1) * (numX + 1)];
            auto& node3 = mesh.Nodes[iX + (iY + 1) * (numX + 1)];
            mesh.Elements.Add({{{node0, node1, node2, node3}, interpolation}});
        }
    return mesh;
}

MeshFem UnitMeshFem::Transform(MeshFem&& oldMesh, std::function<Eigen::VectorXd(Eigen::VectorXd)> f)
{
    // Build a group (MeshFem::NodesTotal() selects all coordinate nodes) to avoid duplicates. Otherwise, the
    // transformation is applied multiple times. This is, however, a bit of a bottleneck for big meshes.
    for (auto& node : oldMesh.NodesTotal())
        node.SetValues(f(node.GetValues()));

    return std::move(oldMesh);
}
