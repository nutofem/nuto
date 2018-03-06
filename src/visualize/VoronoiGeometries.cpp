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

std::vector<double> CenteredPointsWithBoundary(std::vector<double> pts)
{
    std::vector<double> result;
    result.push_back(-1.);
    for (size_t i = 0; i < pts.size() - 1; ++i)
        result.push_back((pts[i] + pts[i + 1]) / 2);
    result.push_back(1.);
    return result;
}

std::vector<double> Gaussspaced(size_t num)
{
    std::vector<double> pts = NuTo::Math::ComputeWeightsAndPoints1DGauss(num).second;
    return CenteredPointsWithBoundary(pts);
}

std::vector<double> Lobattospaced(size_t num)
{
    std::vector<double> pts = NuTo::Math::ComputeWeightsAndPoints1DLobatto(num).second;
    return CenteredPointsWithBoundary(pts);
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

/* =========================================================
 *              TRIANGLE STUFF
 * =========================================================
 */

#include <boost/polygon/voronoi.hpp>
#include <Eigen/Dense> // for Eigen::inverse()

typedef boost::polygon::point_data<double> VoronoiPoint;

int AddUnique(std::vector<Eigen::VectorXd>* rExisingPoints, Eigen::Vector2d newPoint)
{
    for (size_t i = 0; i < rExisingPoints->size(); ++i)
        if ((newPoint - (*rExisingPoints)[i]).squaredNorm() < 1.e-4)
            return i;
    rExisingPoints->push_back(newPoint);
    return rExisingPoints->size() - 1;
}

std::vector<Eigen::Vector2d> TransformedHull(std::vector<Eigen::Vector2d> hullpoints, Eigen::Matrix2d T)
{
    std::vector<Eigen::Vector2d> hull;
    for (Eigen::Vector2d p : hullpoints)
        hull.push_back(T * p);
    return hull;
}

std::vector<VoronoiPoint> TransformedIps(const NuTo::IntegrationTypeBase& integrationType, Eigen::Matrix2d T)
{
    std::vector<VoronoiPoint> points;
    for (int i = 0; i < integrationType.GetNumIntegrationPoints(); ++i)
    {
        Eigen::Vector2d ip = T * Eigen::Vector2d(integrationType.GetLocalIntegrationPointCoordinates(i).x(),
                                                 integrationType.GetLocalIntegrationPointCoordinates(i).y());
        points.push_back({ip.x(), ip.y()});
    }
    return points;
}

std::vector<VoronoiPoint> ReflectPointsOnHull(std::vector<VoronoiPoint> points, std::vector<Eigen::Vector2d> hull)
{
    auto voronoipoints = points;
    for (size_t iEdge = 0; iEdge < hull.size(); ++iEdge)
        for (auto point : points)
        {
            Eigen::Vector2d a = hull[iEdge];
            Eigen::Vector2d b = hull[(iEdge + 1) % hull.size()];
            Eigen::Vector2d n = (b - a).normalized();
            Eigen::Vector2d p(point.x(), point.y());

            Eigen::Vector2d pNew = a + 2 * n * n.dot(p - a) - (p - a);
            voronoipoints.push_back({pNew.x(), pNew.y()});
        }
    return voronoipoints;
}

/* Voronoi diagrams may have infinite edges like:
 *        |
 *     x / \ x
 *      / x \
 *
 * Thus, we have to somehow consider the hull. We simply do that by reflecting all points on that hull
 * [ReflectPointsOnHull] and do voronoi on the original _and_ the reflected points.
 *
 * [Performace is not needed here, since the geometry is calculated once and stored.]
 *
 * x   x  / \  x   x
 *       / x \
 *   x  /     \ x
 *     / x   x \
 *    /_________\
 *
 *       x   x
 *
 *         x
 *
 * This will ensure that the voronoi diagramm has an edge on the hull.
 *
 * boost::polygon::voronoi now creates _cells_. We only evaluate cells that are associated with the origninal points and
 * transform them to a NuTo::Visualize::VoronoiGeometry. [AddUnique] will make sure that we add points at the same
 * coordinate only once.
 */
NuTo::Visualize::VoronoiGeometry TransformedVoronoi(const NuTo::IntegrationTypeBase& integrationType,
                                                    std::vector<Eigen::Vector2d> hullpoints, Eigen::Matrix2d T)
{
    T *= 1.e6; // scaling with a big F to deal with the fact that boost::polygon::voronoi handles ints
    auto hull = TransformedHull(hullpoints, T);
    auto points = TransformedIps(integrationType, T);
    auto voronoipoints = ReflectPointsOnHull(points, hull);

    NuTo::Visualize::VoronoiGeometry geometry;
    // initialize with empty ids and the POLYGON type
    geometry.voronoiCells.resize(points.size(), {{}, NuTo::eCellTypes::POLYGON});

    boost::polygon::voronoi_diagram<double> vd;
    boost::polygon::construct_voronoi(voronoipoints.begin(), voronoipoints.end(), &vd);
    for (const auto& cell : vd.cells())
    {
        if (cell.source_index() >= points.size())
            continue;

        const auto* edge = cell.incident_edge();
        // this edge do while thing is the recommendend way to loop over the edges
        do
        {
            auto v = *edge->vertex0();
            Eigen::Vector2d p = T.inverse() * Eigen::Vector2d(v.x(), v.y()); // transform back

            int pointId = AddUnique(&geometry.pointCoordinates, p);
            geometry.voronoiCells[cell.source_index()].cellCornerIds.push_back(pointId);

            edge = edge->next();
        } while (edge != cell.incident_edge());
    }
    return geometry;
}

VoronoiGeometry NuTo::Visualize::VoronoiGeometryTriangle(const NuTo::IntegrationTypeBase& integrationType)
{
    for (int i = 0; i < integrationType.GetNumIntegrationPoints(); ++i)
        if (integrationType.GetLocalIntegrationPointCoordinates(i).sum() > 1.)
            throw Exception(__PRETTY_FUNCTION__, "Only 2D triangle types are supported");

    Eigen::Vector2d p0(0, 0);
    Eigen::Vector2d p1(1, 0);
    Eigen::Vector2d p2(0, 1);

    // The equilateral system will make the voronoi cells join smoothly with - possibly rotated by 120° - neighbors
    Eigen::Matrix2d ToEquilateral = (Eigen::Matrix2d() << 1, .5, 0, std::sqrt(0.75)).finished();
    return TransformedVoronoi(integrationType, {p0, p1, p2}, ToEquilateral);
}
