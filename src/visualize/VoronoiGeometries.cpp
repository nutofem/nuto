#include "VoronoiGeometries.h"
#include "math/Quadrature.h"

using namespace NuTo::Visualize;

// Helper functions (cpp only)
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

std::vector<double> Linspace(double s, double e, size_t num)
{
    std::vector<double> linspace(num);
    for (size_t i = 0; i < num; ++i)
        linspace[i] = s + (e - s) * i / (num - 1);
    return linspace;
}

std::vector<double> Gaussspaced(size_t num)
{
    std::vector<double> pts = NuTo::Math::ComputeWeightsAndPoints1DGauss(num).second;
    std::vector<double> result;
    result.push_back(-1.);
    for (size_t i = 0; i < pts.size() - 1; ++i)
        result.push_back((pts[i] + pts[i + 1]) / 2);
    result.push_back(1.);
    return result;
}

std::vector<double> Lobattospaced(size_t num)
{
    std::vector<double> pts = NuTo::Math::ComputeWeightsAndPoints1DLobatto(num).second;
    std::vector<double> result;
    result.push_back(-1.);
    for (size_t i = 0; i < pts.size() - 1; ++i)
        result.push_back((pts[i] + pts[i + 1]) / 2);
    result.push_back(1.);
    return result;
}


VoronoiGeometry VoronoiGeometryLineHelper(std::vector<double> gridPoints)
{
    VoronoiGeometry geometry;
    int num = gridPoints.size() - 1;

    for (auto x : gridPoints)
        geometry.pointCoordinates.push_back(Eigen::VectorXd::Constant(1, x));

    for (int i = 0; i < num; i++)
        geometry.voronoiCells.push_back(VoronoiCell{{i, i + 1}, NuTo::eCellTypes::LINE});

    return geometry;
}

VoronoiGeometry VoronoiGeometryQuadHelper(std::vector<double> gridPoints)
{
    VoronoiGeometry geometry;
    int num = gridPoints.size() - 1;

    for (auto y : gridPoints)
        for (auto x : gridPoints)
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

VoronoiGeometry VoronoiGeometryBrickHelper(std::vector<double> gridPoints)
{
    VoronoiGeometry geometry;
    int num = gridPoints.size() - 1;

    for (auto z : gridPoints)
        for (auto y : gridPoints)
            for (auto x : gridPoints)
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

// Member functions
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

VoronoiGeometry NuTo::Visualize::VoronoiGeometryLine(int num, Spacing s)
{
    switch (s)
    {
    case EQUIDISTANT:
        return VoronoiGeometryLineHelper(Linspace(-1, 1, num + 1));
    case LOBATTO:
        return VoronoiGeometryLineHelper(Lobattospaced(num));
    case GAUSS:
        return VoronoiGeometryLineHelper(Gaussspaced(num));
    default:
        throw Exception("Unknown spacing");
    }
}

VoronoiGeometry NuTo::Visualize::VoronoiGeometryQuad(int num, Spacing s)
{
    switch (s)
    {
    case EQUIDISTANT:
        return VoronoiGeometryQuadHelper(Linspace(-1, 1, num + 1));
    case LOBATTO:
        return VoronoiGeometryQuadHelper(Lobattospaced(num));
    case GAUSS:
        return VoronoiGeometryQuadHelper(Gaussspaced(num));
    default:
        throw Exception("Unknown spacing");
    }
}

VoronoiGeometry NuTo::Visualize::VoronoiGeometryBrick(int num, Spacing s)
{
    switch (s)
    {
    case EQUIDISTANT:
        return VoronoiGeometryBrickHelper(Linspace(-1, 1, num + 1));
    case LOBATTO:
        return VoronoiGeometryBrickHelper(Lobattospaced(num));
    case GAUSS:
        return VoronoiGeometryBrickHelper(Gaussspaced(num));
    default:
        throw Exception("Unknown spacing");
    }
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
