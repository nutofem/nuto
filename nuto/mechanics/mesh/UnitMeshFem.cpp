#include "nuto/mechanics/mesh/UnitMeshFem.h"
#include "nuto/mechanics/interpolation/InterpolationTriangleLinear.h"
#include "nuto/mechanics/interpolation/InterpolationTrussLinear.h"
#include "nuto/mechanics/interpolation/InterpolationQuadLinear.h"
#include "nuto/mechanics/interpolation/InterpolationBrickLinear.h"

using namespace NuTo;

// anonymous helper functions for mesh creation and transformation
MeshFem CreateNodes1D(int numX)
{
    MeshFem mesh;
    for (int iX = 0; iX < numX + 1; ++iX)
    {
        const double x = static_cast<double>(iX) / numX;
        mesh.CoordinateNodes.Add({x});
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
            mesh.CoordinateNodes.Add(Eigen::Vector2d(x, y));
        }
    return mesh;
}


MeshFem UnitMeshFem::CreateLines(int numX)
{
    MeshFem mesh = CreateNodes1D(numX);
    const auto& interpolation = mesh.CreateInterpolation(NuTo::InterpolationTrussLinear());
    for (int i = 0; i < numX; ++i)
    {
        auto& nl = mesh.CoordinateNodes[i];
        auto& nr = mesh.CoordinateNodes[i + 1];
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
            auto& node0 = mesh.CoordinateNodes[iX + iY * (numX + 1)];
            auto& node1 = mesh.CoordinateNodes[iX + 1 + iY * (numX + 1)];
            auto& node2 = mesh.CoordinateNodes[iX + 1 + (iY + 1) * (numX + 1)];
            auto& node3 = mesh.CoordinateNodes[iX + (iY + 1) * (numX + 1)];
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
            auto& node0 = mesh.CoordinateNodes[iX + iY * (numX + 1)];
            auto& node1 = mesh.CoordinateNodes[iX + 1 + iY * (numX + 1)];
            auto& node2 = mesh.CoordinateNodes[iX + 1 + (iY + 1) * (numX + 1)];
            auto& node3 = mesh.CoordinateNodes[iX + (iY + 1) * (numX + 1)];
            mesh.Elements.Add({{{node0, node1, node2, node3}, interpolation}});
        }
    return mesh;
}

MeshFem UnitMeshFem::CreateBricks(int numX, int numY, int numZ)
{
    MeshFem mesh;
    int numXe = numX + 1;
    int numYe = numY + 1;
    int numZe = numZ + 1;
    for (int iZ = 0; iZ < numZe; ++iZ)
        for (int iY = 0; iY < numYe; ++iY)
            for (int iX = 0; iX < numXe; ++iX)
            {
                const double x = static_cast<double>(iX) / numX;
                const double y = static_cast<double>(iY) / numY;
                const double z = static_cast<double>(iZ) / numZ;
                mesh.CoordinateNodes.Add(Eigen::Vector3d(x, y, z));
            }
    const auto& interpolation = mesh.CreateInterpolation(NuTo::InterpolationBrickLinear());
    for (int iZ = 0; iZ < numZ; ++iZ)
        for (int iY = 0; iY < numY; ++iY)
            for (int iX = 0; iX < numX; ++iX)
            {
                auto& node0 = mesh.CoordinateNodes[iX + iY * numXe + iZ * numXe * numYe];
                auto& node1 = mesh.CoordinateNodes[iX + 1 + iY * numXe + iZ * numXe * numYe];
                auto& node2 = mesh.CoordinateNodes[iX + 1 + (iY + 1) * numXe + iZ * numXe * numYe];
                auto& node3 = mesh.CoordinateNodes[iX + (iY + 1) * numXe + iZ * numXe * numYe];
                auto& node4 = mesh.CoordinateNodes[iX + iY * numXe + (iZ + 1) * numXe * numYe];
                auto& node5 = mesh.CoordinateNodes[iX + 1 + iY * numXe + (iZ + 1) * numXe * numYe];
                auto& node6 = mesh.CoordinateNodes[iX + 1 + (iY + 1) * numXe + (iZ + 1) * numXe * numYe];
                auto& node7 = mesh.CoordinateNodes[iX + (iY + 1) * numXe + (iZ + 1) * numXe * numYe];
                mesh.Elements.Add({{{node0, node1, node2, node3, node4, node5, node6, node7}, interpolation}});
            }
    return mesh;
}


MeshFem UnitMeshFem::Transform(MeshFem&& oldMesh, std::function<Eigen::VectorXd(Eigen::VectorXd)> f)
{
    // Build a group (MeshFem::NodesTotal() selects all coordinate nodes) to avoid duplicates. Otherwise, the
    // transformation is applied multiple times. This is, however, a bit of a bottleneck for big meshes.
    for (auto& node : oldMesh.NodesTotal())
        node.SetCoordinates(f(node.GetCoordinates()));

    return std::move(oldMesh);
}
