/*
 * DelaunayVoronoi.h
 *
 *  Created on: 1 Sep 2015
 *      Author: ttitsche
 */
#pragma once

#include <iostream>
#include <Eigen/Dense>

#include <array>
#include <list>
#include <set>
#include <map>
#include <vector>
#include <unordered_set>
#include <unordered_map>
#include <algorithm>
#include <cmath>

#include <string>

#include <fstream>

class Edge
{
public:
    Edge()
    {
    }

    Edge(Eigen::Vector2d rP0, Eigen::Vector2d rP1)
    {
        mPoints[0] = rP0;
        mPoints[1] = rP1;
    }

    const Eigen::Vector2d P0() const
    {
        return mPoints[0];
    }

    const Eigen::Vector2d P1() const
    {
        return mPoints[1];
    }

    bool operator==(const Edge& rOther) const
    {
        bool result = (P0() == rOther.P0() && P1() == rOther.P1()) || (P0() == rOther.P1() && P1() == rOther.P0());
        return result;
    }

    bool Contains(Eigen::Vector2d rPoint) const
    {
        if (P0() == rPoint || P1() == rPoint)
            return true;
        return false;
    }

    void Info() const
    {
        std::cout << "E: (" << P0().transpose() << ") and (" << P1().transpose() << ") \n";
    }

    size_t GetHash() const
    {
        return std::hash<double>()(P0().x()) + std::hash<double>()(P0().y()) + std::hash<double>()(P1().x()) +
               std::hash<double>()(P1().y());
    }

    //! @brief returns the line parameter l of the line rPoint + l * rDirection
    //! where the line intersects with the edge. Returns -1 if there is no intersection.
    double GetIntersectionLineParameter(const Eigen::Vector2d rPoint, const Eigen::Vector2d rDirection) const
    {
        const Eigen::Vector2d e = mPoints[1] - mPoints[0];
        const Eigen::Vector2d& n = rDirection;

        double eXn = e.x() * n.y() - e.y() * n.x(); // pseudo cross product 2D
        double nXe = -eXn;

        const Eigen::Vector2d D = rPoint - mPoints[0];

        double edgeParameter = (D.x() * n.y() - D.y() * n.x()) / eXn;
        double lineParameter = (D.y() * e.x() - D.x() * e.y()) / nXe;

        if (edgeParameter < 0. || edgeParameter > 1.)
            return -1.;

        return lineParameter;
    }

    //! @brief calculates the (perpendicular) distance of rPoint to the edge
    double GetDistance(Eigen::Vector2d rPoint) const
    {
        Eigen::Vector2d p = P1() - P0();
        Eigen::Vector2d a = rPoint - P0();

        return (a - a.dot(p) * p / p.dot(p)).norm();
    }

private:
    std::array<Eigen::Vector2d, 2> mPoints;
};

class Triangle
{
public:
    Triangle(Edge rE0, Edge rE1, Edge rE2)
    {
        mEdges[0] = rE0;
        mEdges[1] = rE1;
        mEdges[2] = rE2;
    }

    std::array<Eigen::Vector2d, 3> GetUniquePoints() const
    {
        std::array<Eigen::Vector2d, 3> uniquePoints;
        uniquePoints[0] = mEdges[0].P0();
        uniquePoints[1] = mEdges[0].P1();

        if (mEdges[0].Contains(mEdges[1].P0()))
            uniquePoints[2] = mEdges[1].P1();
        else
            uniquePoints[2] = mEdges[1].P0();

        assert(uniquePoints[0] != uniquePoints[1]);
        assert(uniquePoints[0] != uniquePoints[2]);
        assert(uniquePoints[1] != uniquePoints[2]);

        return uniquePoints;
    }

    std::array<Edge, 3> GetEdges() const
    {
        return mEdges;
    }

    bool Contains(const Eigen::Vector2d& rPoint) const
    {
        for (auto edge : mEdges)
            if (edge.Contains(rPoint))
                return true;
        return false;
    }

    bool Contains(Edge rEdge) const
    {
        for (auto edge : mEdges)
            if (rEdge == edge)
                return true;
        return false;
    }

    bool operator==(const Triangle& rOther) const
    {
        for (auto otherPoint : rOther.GetUniquePoints())
            if (not Contains(otherPoint))
                return false;

        return true;
    }

    // https://en.wikipedia.org/wiki/Circumscribed_circle#Cartesian_coordinates_2
    Eigen::Vector2d GetCircumcenter() const
    {
        auto points = GetUniquePoints();
        Eigen::Vector2d A = points[0];
        Eigen::Vector2d B = points[1];
        Eigen::Vector2d C = points[2];

        double d = 2 * (A.x() * (B.y() - C.y()) + B.x() * (C.y() - A.y()) + C.x() * (A.y() - B.y()));

        double x = A.dot(A) * (B.y() - C.y()) + B.dot(B) * (C.y() - A.y()) + C.dot(C) * (A.y() - B.y());
        double y = A.dot(A) * (C.x() - B.x()) + B.dot(B) * (A.x() - C.x()) + C.dot(C) * (B.x() - A.x());

        return Eigen::Vector2d(x / d, y / d);
    }

    bool CircumcircleContainsPoint(const Eigen::Vector2d& rPoint) const
    {
        auto points = GetUniquePoints();
        Eigen::Vector2d A = points[0];
        Eigen::Vector2d B = points[1];
        Eigen::Vector2d C = points[2];
        Eigen::Vector2d D = rPoint;

        // check if clockwise oriented A --> B --> C ...
        Eigen::Vector2d v1 = B - A;
        Eigen::Vector2d v2 = C - A;

        double dir = v1.x() * v2.y() - v1.y() * v2.x();
        // ... and swap B and C if not.
        if (dir < 0.)
        {
            auto q = B;
            B = C;
            C = q;
        }

        if (A == D || B == D || C == D)
            return false;

        Eigen::Matrix3d matr;

        matr(0, 0) = A.x() - D.x();
        matr(1, 0) = B.x() - D.x();
        matr(2, 0) = C.x() - D.x();

        matr(0, 1) = A.y() - D.y();
        matr(1, 1) = B.y() - D.y();
        matr(2, 1) = C.y() - D.y();

        matr(0, 2) = A.x() * A.x() - D.x() * D.x() + A.y() * A.y() - D.y() * D.y();
        matr(1, 2) = B.x() * B.x() - D.x() * D.x() + B.y() * B.y() - D.y() * D.y();
        matr(2, 2) = C.x() * C.x() - D.x() * D.x() + C.y() * C.y() - D.y() * D.y();

        return matr.determinant() > 0;
    }

    void Info() const
    {
        auto points = GetUniquePoints();
        std::cout << "T: \n"
                  << points[0].transpose() << " | " << points[1].transpose() << " | " << points[2].transpose()
                  << " | \n";
    }

private:
    std::array<Edge, 3> mEdges;
};

class Polygon
{
public:
    Polygon()
    {
    }

    Polygon(std::vector<Eigen::Vector2d> rPoints)
        : mPoints(rPoints)
    {
        SortAndRemoveDoublicates();
    }

    std::vector<Eigen::Vector2d> GetPoints()
    {
        return mPoints;
    }

    void SortAndRemoveDoublicates()
    {
        // find some point inside the (convex) voronoi polygon
        Eigen::Vector2d center(0., 0.);
        for (auto point : mPoints)
            center += point;
        center /= mPoints.size();

        std::map<double, Eigen::Vector2d> anglesPoints;

        // get angle of each voronoi point
        for (auto vp : mPoints)
        {
            double angle = atan2(vp.x() - center.x(), vp.y() - center.y());
            anglesPoints[angle] = vp;
        }

        std::vector<Eigen::Vector2d> sortedUniqueVoronoiPoints;

        // two points with an angle difference < delta
        // are the same points and are only considered once.
        double delta = 1.e-6;
        double lastAngle = -10;
        for (auto pair : anglesPoints)
        {
            double currentAngle = pair.first;
            if (currentAngle - lastAngle > delta)
            {
                sortedUniqueVoronoiPoints.push_back(pair.second);
                lastAngle = currentAngle;
            }
        }
        mPoints = sortedUniqueVoronoiPoints;
        UpdateEdges();
    }

    void UpdateEdges()
    {
        mEdges.clear();
        int numPoints = mPoints.size();
        for (int i = 0; i < numPoints; ++i)
            mEdges.push_back(Edge(mPoints[i], mPoints[(i + 1) % numPoints]));
    }

    //! @brief counts the intersections with the edges. The point rPoint is inside the polygon if that number is odd.
    bool PointIsInside(Eigen::Vector2d rPoint) const
    {
        int intersectionCounter = 0;
        Eigen::Vector2d anyDirection(1., 0.);

        for (auto edge : mEdges)
            if (edge.GetIntersectionLineParameter(rPoint, anyDirection) > 0.)
                intersectionCounter++;

        return intersectionCounter % 2 == 1;
    }

    //! @brief returns true, if the polygon intersects the segment rP0-rP1
    bool SegmentIntersects(Eigen::Vector2d rP0, Eigen::Vector2d rP1) const
    {
        for (auto edge : mEdges)
        {
            double lineParameter = edge.GetIntersectionLineParameter(rP0, rP1 - rP0);
            if (lineParameter > 0. and lineParameter < 1.)
                return true;
        }
        return false;
    }

    //! @brief A line from rPoint in rDirection may intersect multiple times with the polygon. This method returns the
    //! closest one to rPoint.
    Eigen::Vector2d GetClosestIntersection(Eigen::Vector2d rPoint, Eigen::Vector2d rDirection) const
    {
        std::vector<Eigen::Vector2d> intersections;
        double closestLineParameter = 1.e10;

        for (auto edge : mEdges)
        {
            double lineParameter = edge.GetIntersectionLineParameter(rPoint, rDirection);
            if (lineParameter > 1.e-10)
                closestLineParameter = std::min(closestLineParameter, lineParameter);
        }
        //        if (closestLineParameter == 1.e10)
        //            throw "No intersection found!";

        return Eigen::Vector2d(rPoint + rDirection * closestLineParameter);
    }

    //! @brief Returns the index of the edge that contains rPoint
    int GetEdgeIndex(Eigen::Vector2d rPoint) const
    {
        for (unsigned int i = 0; i < mEdges.size(); ++i)
            if (mEdges[i].GetDistance(rPoint) < 1.e-8)
                return i;
        return -1;
    }

    //! @brief Returns a boundary point that is contained by both rEdge1 and rEdge2.
    Eigen::Vector2d GetPointOnBothEdges(int rEdge1, int rEdge2) const
    {
        auto e1 = mEdges[rEdge1];
        auto e2 = mEdges[rEdge2];
        if (e2.Contains(e1.P0()))
            return e1.P0();
        if (e2.Contains(e1.P1()))
            return e1.P1();

        throw "There is no point on both rEdge1 and rEdge2";
    }

    //! @brief returns all mPoints except rPoint
    std::vector<Eigen::Vector2d> GetOtherPointsOnBoundary(Eigen::Vector2d rPoint) const
    {
        std::vector<Eigen::Vector2d> points = mPoints;
        auto it = std::find(points.begin(), points.end(), rPoint);
        assert(it != points.end());

        points.erase(it);
        return points;
    }

    //! @brief add points to this polygon if it intersects with the boundary
    //! If the polygon contains two neighboring points that are located on different edges of rBoundary
    //! the polygon must include at least one point of rBoundary. Read further comments...
    void AddIntersectionWithBoundary(const Polygon& rBoundary, const Eigen::Vector2d& rIP)
    {

        int numPointsOnBoundary = 0;
        int numPoints = mPoints.size();

        for (auto vp : mPoints)
            if (rBoundary.GetEdgeIndex(vp) != -1)
                numPointsOnBoundary++;

        if (numPointsOnBoundary != 2)
            return;


        // may contain a boundary point
        // needed if two neigboring points are on different boundary edges
        for (int i = 0; i < numPoints; ++i)
        {
            int e0 = rBoundary.GetEdgeIndex(mPoints[i]);
            int e1 = rBoundary.GetEdgeIndex(mPoints[(i + 1) % numPoints]);

            if (e0 == -1 || e1 == -1 || e0 == e1)
                continue;

            auto boundaryPoint = rBoundary.GetPointOnBothEdges(e0, e1);

            // the boundaryPoint must be inside the voronoi cell
            auto originalPoints = mPoints; // backup

            mPoints.push_back(boundaryPoint);
            SortAndRemoveDoublicates();

            if (PointIsInside(rIP))
            {
                break; // since only one edge per voronoi polygon allowed
                // normal case
            }
            else
            {
                std::cout << "This happens!" << std::endl;
                // restore points
                mPoints = originalPoints;

                auto otherBoundaryPoints = rBoundary.GetOtherPointsOnBoundary(boundaryPoint);
                for (auto point : otherBoundaryPoints)
                    mPoints.push_back(point);

                SortAndRemoveDoublicates();
                //                if (not PointIsInside(rIP))
                //                    throw "Damn it!";
            }
        }
    }

private:
    std::vector<Edge> mEdges;
    std::vector<Eigen::Vector2d> mPoints;
};


class DelaunayVoronoi
{
public:
    DelaunayVoronoi(const std::vector<Eigen::Vector2d>& rPoints, bool rCalculateInTransformedSystem)
        : mPoints(rPoints)
        , mCalculateInTransformedSystem(rCalculateInTransformedSystem)
    {
        if (mCalculateInTransformedSystem)
            Transform(mPoints);
    }

    void SetBoundary(const std::vector<Eigen::Vector2d>& rBoundaryPoints)
    {
        auto points = rBoundaryPoints;
        if (mCalculateInTransformedSystem)
            Transform(points);

        mBoundary = Polygon(points);
    }

    std::list<Triangle> GetDelaunayTriangulation()
    {
        if (mTriangulation.size() == 0)
            CalculateDelaunayTriangulation();

        return mTriangulation;
    }

    std::vector<Polygon> GetVoronoiCells()
    {
        if (mVoronoiPolygons.size() == 0)
            CalculateVoronoiPolygons();

        return mVoronoiPolygons;
    }

    void Transform(std::vector<Eigen::Vector2d>& rPoints)
    {
        Eigen::Matrix2d T;
        T << 1, .5, 0, std::sqrt(0.75);

        for (auto& point : rPoints)
        {
            auto p = point;
            point = T * p;
        }
    }

    void TransformInverse(std::vector<Eigen::Vector2d>& rPoints)
    {
        Eigen::Matrix2d T;
        T << 1, .5, 0, std::sqrt(0.75);

        auto invT = T.inverse();

        for (auto& point : rPoints)
        {
            auto p = point;
            point = invT * p;
        }
    }

    //! @brief special thanks to https://en.wikipedia.org/wiki/Bowyer%E2%80%93Watson_algorithm
    void CalculateDelaunayTriangulation()
    {
        mTriangulation.clear();

        // define a super quad as two super triangles that contain all points and an offset
        Eigen::Vector2d superMin = mPoints[0];
        Eigen::Vector2d superMax = mPoints[0];
        for (auto point : mPoints)
        {
            superMin = superMin.cwiseMin(point);
            superMax = superMax.cwiseMax(point);
        }

        auto size = superMax - superMin;
        superMin -= 3 * size;
        superMax += 3 * size;

        std::array<Eigen::Vector2d, 4> sP;
        sP[0] = superMin;
        sP[1] = Eigen::Vector2d(superMax.x(), superMin.y());
        sP[2] = superMax;
        sP[3] = Eigen::Vector2d(superMin.x(), superMax.y());

        mTriangulation.push_back(Triangle(Edge(sP[0], sP[1]), Edge(sP[1], sP[2]), Edge(sP[2], sP[0])));
        mTriangulation.push_back(Triangle(Edge(sP[0], sP[2]), Edge(sP[2], sP[3]), Edge(sP[3], sP[0])));


        // add every point from rPoints using the Bowyer–Watson algorithm
        for (auto point : mPoints)
        {
            std::list<const Triangle*> badTriangles;
            for (Triangle& triangle : mTriangulation)
            {
                if (triangle.CircumcircleContainsPoint(point))
                    badTriangles.push_back(&triangle);
            }

            std::list<Edge> polygon;
            for (const Triangle* badTriangle : badTriangles)
            {
                for (Edge edge : badTriangle->GetEdges())
                {
                    bool edgeNotShared = true;
                    for (const Triangle* otherBadTriangle : badTriangles)
                    {
                        if (*otherBadTriangle == *badTriangle)
                            continue;

                        if (otherBadTriangle->Contains(edge))
                        {
                            edgeNotShared = false;
                            break;
                        }
                    }
                    if (edgeNotShared)
                        polygon.push_back(edge);
                }
            }

            for (const Triangle* triangle : badTriangles)
                mTriangulation.remove(*triangle);

            // create new edges containing the point
            for (Edge edge : polygon)
                mTriangulation.push_back(Triangle(edge, Edge(point, edge.P0()), Edge(point, edge.P1())));
        }

        std::list<const Triangle*> trianglesToErase;

        // remove triangles containing the super points
        for (Triangle& triangle : mTriangulation)
            for (auto superPoint : sP)
                if (triangle.Contains(superPoint))
                    trianglesToErase.push_back(&triangle);

        for (const Triangle* triangleToErase : trianglesToErase)
            mTriangulation.remove(*triangleToErase);


        std::ofstream file;
        file.open("triangles.dat");
        for (auto triangle : mTriangulation)
        {
            auto points = triangle.GetUniquePoints();
            file << points[0].transpose() << '\n'
                 << points[1].transpose() << '\n'
                 << points[2].transpose() << '\n'
                 << points[0].transpose() << "\n \n";
        }
        file.close();
    }


    void CalculateVoronoiPolygons()
    {
        if (mBoundary.GetPoints().size() == 0)
            throw "Define a the boundary first! [AddBoundary()]!";

        int num = mPoints.size();
        mVoronoiPolygons.resize(num);

        if (mTriangulation.size() == 0)
            CalculateDelaunayTriangulation();


        struct EdgeHash
        {
            size_t operator()(const Edge& x) const
            {
                return x.GetHash();
            }
        };

        // create a data structure that contains all the triangles for every edge
        std::unordered_map<Edge, std::vector<const Triangle*>, EdgeHash> edgeToTriangle;
        for (const Triangle& triangle : mTriangulation)
            for (Edge edge : triangle.GetEdges())
                edgeToTriangle[edge].push_back(&triangle);

        for (int p = 0; p < num; ++p)
        {
            auto point = mPoints[p];

            std::list<Edge> edges;
            for (auto pair : edgeToTriangle)
            {
                const Edge& edge = pair.first;
                if (edge.Contains(point))
                    edges.push_back(edge);
            }

            std::vector<Eigen::Vector2d> voronoiPoints;

            for (const Edge& edge : edges)
            {
                std::vector<const Triangle*> triangles = edgeToTriangle[edge];
                assert(triangles.size() > 0);
                assert(triangles.size() <= 2);

                Eigen::Vector2d p0, p1;

                if (triangles.size() == 1)
                {

                    p0 = triangles[0]->GetCircumcenter();

                    // build a vector n from the circumcenter to its edge ...
                    Eigen::Vector2d e = edge.P1() - edge.P0();
                    Eigen::Vector2d edgePoint = edge.P0() + e.dot(p0 - edge.P0()) / e.dot(e) * e;
                    Eigen::Vector2d n = edgePoint - p0;

                    // if the circumcenter lays on the edge, n is 0.
                    // instead pick any n perpendicular to e;

                    if (std::abs(n.dot(n)) < 1.e-8)
                        n = Eigen::Vector2d(e[1], -e[0]);

                    // For the special case of a circumcenter that is
                    // outside of its triangle and inside the boundary:
                    //
                    // a line l = p0 + t * n with t > 1 must not intersect
                    // any other edges of the triangle.
                    // if it does, reverse: n --> -n
                    bool intersectsAgain = false;

                    for (Edge edgeTri : triangles[0]->GetEdges())
                    {
                        if (edgeTri == edge)
                            continue;

                        double t = edgeTri.GetIntersectionLineParameter(p0, n);
                        if (t > 1.)
                            intersectsAgain = true;
                    }

                    if (intersectsAgain)
                        n = -n;

                    if (not mBoundary.PointIsInside(p0))
                    {

                        if (intersectsAgain)
                            // do not consider the edge
                            // that intersects the triangle again
                            // oisch... sketches would help here... believe me please :)
                            continue;

                        p0 = mBoundary.GetClosestIntersection(p0, n);
                    }
                    p1 = mBoundary.GetClosestIntersection(p0, n);
                }

                if (triangles.size() == 2)
                {
                    // find the point p0 in the boundary
                    p0 = triangles[0]->GetCircumcenter();
                    p1 = triangles[1]->GetCircumcenter();

                    if (not mBoundary.PointIsInside(p0))
                    {
                        // if the line from p0 to p1 does not intersect
                        // the boundary, continue
                        if (not mBoundary.SegmentIntersects(p0, p1))
                            continue;

                        // find p0 as closest intersection of p0 to p1
                        p0 = mBoundary.GetClosestIntersection(p0, p1 - p0);
                    }

                    if (not mBoundary.PointIsInside(p1))
                    {
                        // find the intersection of the line p0-p1 with
                        // the boundary
                        p1 = mBoundary.GetClosestIntersection(p1, p0 - p1);
                    }
                }

                voronoiPoints.push_back(p0);
                voronoiPoints.push_back(p1);
            }

            Polygon voronoiCell(voronoiPoints);
            voronoiCell.AddIntersectionWithBoundary(mBoundary, point);
            mVoronoiPolygons[p] = voronoiCell;
        }
    }

    void CalculateVisualizationCellsPolygon(std::vector<Eigen::Vector2d>& rVisuPoints,
                                            std::vector<std::vector<unsigned int>>& rVisuCellIndices)
    {
        if (mVoronoiPolygons.size() == 0)
            CalculateVoronoiPolygons();

        rVisuPoints.clear();
        rVisuCellIndices.resize(mPoints.size());

        for (auto polygon : mVoronoiPolygons)
            for (auto point : polygon.GetPoints())
                if (std::find(rVisuPoints.begin(), rVisuPoints.end(), point) == rVisuPoints.end())
                    rVisuPoints.push_back(point);

        for (unsigned int iIP = 0; iIP < mPoints.size(); ++iIP)
        {
            auto polygon = mVoronoiPolygons[iIP];

            // perform a triangulation of the convex polygon
            std::vector<unsigned int> visuCell;
            for (auto vPoint : polygon.GetPoints())
                // subtracting two iterators returns the index
                visuCell.push_back(std::find(rVisuPoints.begin(), rVisuPoints.end(), vPoint) - rVisuPoints.begin());
            rVisuCellIndices[iIP] = visuCell;
        }

        if (mCalculateInTransformedSystem)
            TransformInverse(rVisuPoints);
    }

    //! @brief calculates the polygon visualization cells and triangulates them again
    //! @param rVisuPoints ... vector of unique visualization points
    //! @param rVisuCellIndices ... vector of triangular visualization cells. Each cell contains a vector of visuPoint
    //! indices
    //! @param rCellIPIndex ... contains the index of the integration point of each visualization cell
    void CalculateVisualizationCellsTriangle(std::vector<Eigen::Vector2d>& rVisuPoints,
                                             std::vector<std::array<unsigned int, 3>>& rVisuCellIndices,
                                             std::vector<unsigned int>& rCellIPIndex)
    {
        rVisuCellIndices.clear();
        rCellIPIndex.clear();

        std::vector<std::vector<unsigned int>> polygonCells;
        CalculateVisualizationCellsPolygon(rVisuPoints, polygonCells);
        for (unsigned int iIP = 0; iIP < polygonCells.size(); ++iIP)
        {
            auto& polygon = polygonCells[iIP];

            unsigned int p0 = polygon[0];
            // split each polygonCell in multiple triangle cells
            for (unsigned int i = 1; i < polygon.size() - 1; ++i)
            {
                unsigned int p1 = polygon[i];
                unsigned int p2 = polygon[i + 1];

                rVisuCellIndices.push_back({{p0, p1, p2}});
                rCellIPIndex.push_back(iIP);
            }
        }
    }


private:
    std::vector<Eigen::Vector2d> mPoints;
    std::list<Triangle> mTriangulation;
    std::vector<Polygon> mVoronoiPolygons;
    Polygon mBoundary;
    bool mCalculateInTransformedSystem;
};
