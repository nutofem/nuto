#include "VoronoiGeometries.h"

using namespace NuTo::Visualize;

std::vector<double> Linspace(double s, double e, size_t num)
{
    std::vector<double> linspace(num);
    for (size_t i = 0; i < num; ++i)
        linspace[i] = s + (e - s) * i / (num - 1);
    return linspace;
}

VoronoiGeometry NuTo::Visualize::VoronoiGeometryLine(int num)
{
    VoronoiGeometry geometry;
    for (auto x : Linspace(-1, 1, num + 1))
        geometry.pointCoordinates.push_back(Eigen::VectorXd::Constant(1, x));

    for (int i = 0; i < num; i++)
        geometry.voronoiCells.push_back(VoronoiCell{{i, i + 1}, NuTo::eCellTypes::LINE});

    return geometry;
}

VoronoiGeometry NuTo::Visualize::VoronoiGeometryQuad(int num)
{
    VoronoiGeometry geometry;
    for (auto y : Linspace(-1, 1, num + 1))
        for (auto x : Linspace(-1, 1, num + 1))
            geometry.pointCoordinates.push_back(Eigen::Vector2d(x, y));

    for (int row = 0; row < num; row++)
    {
        for (int col = 0; col < num; col++)
        {
            int start = row * (num + 1) + col;
            geometry.voronoiCells.push_back(
                    VoronoiCell{{start, start + 1, start + num + 2, start + num + 1}, NuTo::eCellTypes::QUAD});
        }
    }
    return geometry;
}

VoronoiGeometry NuTo::Visualize::VoronoiGeometryBrick(int num)
{
    VoronoiGeometry geometry;
    for (auto z : Linspace(-1, 1, num + 1))
        for (auto y : Linspace(-1, 1, num + 1))
            for (auto x : Linspace(-1, 1, num + 1))
                geometry.pointCoordinates.push_back(Eigen::Vector3d(x, y, z));

    for (int height = 0; height < num; height++)
    {
        for (int row = 0; row < num; row++)
        {
            for (int col = 0; col < num; col++)
            {
                int start1 = row * (num + 1) + col + height * ((num + 1) * (num + 1));
                int start2 = row * (num + 1) + col + (height + 1) * ((num + 1) * (num + 1));
                geometry.voronoiCells.push_back(VoronoiCell{{start1, start1 + 1, start1 + num + 2, start1 + num + 1,
                                                             start2, start2 + 1, start2 + num + 2, start2 + num + 1},
                                                            NuTo::eCellTypes::HEXAHEDRON});
            }
        }
    }
    return geometry;
}

VoronoiGeometry NuTo::Visualize::VoronoiGeometryTriangle3Ip()
{
    VoronoiGeometry geometry;
    geometry.pointCoordinates.push_back(Eigen::Vector2d(0, 0));
    geometry.pointCoordinates.push_back(Eigen::Vector2d(0.5, 0));
    geometry.pointCoordinates.push_back(Eigen::Vector2d(1, 0));
    geometry.pointCoordinates.push_back(Eigen::Vector2d(0, 0.5));
    geometry.pointCoordinates.push_back(Eigen::Vector2d(1. / 3., 1. / 3.));
    geometry.pointCoordinates.push_back(Eigen::Vector2d(0.5, 0.5));
    geometry.pointCoordinates.push_back(Eigen::Vector2d(0, 1));

    geometry.voronoiCells.push_back({{0, 1, 4, 3}, NuTo::eCellTypes::QUAD});
    geometry.voronoiCells.push_back({{1, 2, 5, 4}, NuTo::eCellTypes::QUAD});
    geometry.voronoiCells.push_back({{4, 5, 6, 3}, NuTo::eCellTypes::QUAD});

    return geometry;
}
