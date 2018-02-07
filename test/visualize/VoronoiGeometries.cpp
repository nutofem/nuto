#include <iostream>
#include "BoostUnitTest.h"
#include "visualize/VoronoiGeometries.h"

using namespace NuTo;

BOOST_AUTO_TEST_CASE(VoronoiGeometryLine)
{
    int numDivisions = 3;

    auto v1 = Visualize::VoronoiGeometryLine(numDivisions, Visualize::EQUIDISTANT);
    auto v2 = Visualize::VoronoiGeometryLine(numDivisions, Visualize::LOBATTO);
    auto v3 = Visualize::VoronoiGeometryLine(numDivisions, Visualize::GAUSS);

    BOOST_CHECK_EQUAL(v1.pointCoordinates.size(), numDivisions + 1);
    BOOST_CHECK_EQUAL(v2.pointCoordinates.size(), numDivisions + 1);
    BOOST_CHECK_EQUAL(v3.pointCoordinates.size(), numDivisions + 1);
}

BOOST_AUTO_TEST_CASE(VoronoiGeometryQuad)
{
    int numDivisions = 3;
    int numPoints = (numDivisions + 1) * (numDivisions + 1);

    auto v1 = Visualize::VoronoiGeometryQuad(numDivisions, Visualize::EQUIDISTANT);
    auto v2 = Visualize::VoronoiGeometryQuad(numDivisions, Visualize::LOBATTO);
    auto v3 = Visualize::VoronoiGeometryQuad(numDivisions, Visualize::GAUSS);

    BOOST_CHECK_EQUAL(v1.pointCoordinates.size(), numPoints);
    BOOST_CHECK_EQUAL(v2.pointCoordinates.size(), numPoints);
    BOOST_CHECK_EQUAL(v3.pointCoordinates.size(), numPoints);
}

BOOST_AUTO_TEST_CASE(VoronoiGeometryBrick)
{
    int numDivisions = 3;
    int numPoints = (numDivisions + 1) * (numDivisions + 1) * (numDivisions + 1);

    auto v1 = Visualize::VoronoiGeometryBrick(numDivisions, Visualize::EQUIDISTANT);
    auto v2 = Visualize::VoronoiGeometryBrick(numDivisions, Visualize::LOBATTO);
    auto v3 = Visualize::VoronoiGeometryBrick(numDivisions, Visualize::GAUSS);

    BOOST_CHECK_EQUAL(v1.pointCoordinates.size(), numPoints);
    BOOST_CHECK_EQUAL(v2.pointCoordinates.size(), numPoints);
    BOOST_CHECK_EQUAL(v3.pointCoordinates.size(), numPoints);
}
