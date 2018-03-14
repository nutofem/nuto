#include <iostream>
#include "BoostUnitTest.h"
#include "visualize/VoronoiGeometries.h"
#include "mechanics/integrationtypes/IntegrationTypeTensorProduct.h"
#include "mechanics/integrationtypes/IntegrationTypeTetrahedron.h"
#include "mechanics/integrationtypes/IntegrationTypeTriangle.h"

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

BOOST_AUTO_TEST_CASE(Voronoi2D_1)
{
    auto v = Visualize::VoronoiGeometryQuad(1);
    BOOST_CHECK_EQUAL(v.pointCoordinates.size(), 4);
    BoostUnitTest::CheckEigenMatrix(v.pointCoordinates[0], Eigen::Vector2d(-1, -1));
    BoostUnitTest::CheckEigenMatrix(v.pointCoordinates[1], Eigen::Vector2d(1, -1));
    BoostUnitTest::CheckEigenMatrix(v.pointCoordinates[2], Eigen::Vector2d(-1, 1));
    BoostUnitTest::CheckEigenMatrix(v.pointCoordinates[3], Eigen::Vector2d(1, 1));

    BoostUnitTest::CheckVector(v.voronoiCells[0].cellCornerIds, std::vector<double>{0, 1, 3, 2}, 4);
}

BOOST_AUTO_TEST_CASE(Voronoi3D_1)
{
    auto v = Visualize::VoronoiGeometryBrick(1);
    BOOST_CHECK_EQUAL(v.pointCoordinates.size(), 8);
    BoostUnitTest::CheckEigenMatrix(v.pointCoordinates[0], Eigen::Vector3d(-1, -1, -1));
    BoostUnitTest::CheckEigenMatrix(v.pointCoordinates[1], Eigen::Vector3d(1, -1, -1));
    BoostUnitTest::CheckEigenMatrix(v.pointCoordinates[2], Eigen::Vector3d(-1, 1, -1));
    BoostUnitTest::CheckEigenMatrix(v.pointCoordinates[3], Eigen::Vector3d(1, 1, -1));
    BoostUnitTest::CheckEigenMatrix(v.pointCoordinates[4], Eigen::Vector3d(-1, -1, 1));
    BoostUnitTest::CheckEigenMatrix(v.pointCoordinates[5], Eigen::Vector3d(1, -1, 1));
    BoostUnitTest::CheckEigenMatrix(v.pointCoordinates[6], Eigen::Vector3d(-1, 1, 1));
    BoostUnitTest::CheckEigenMatrix(v.pointCoordinates[7], Eigen::Vector3d(1, 1, 1));

    auto vCellIds = v.voronoiCells[0].cellCornerIds;
    BOOST_CHECK_EQUAL(vCellIds[0], 0);
    BOOST_CHECK_EQUAL(vCellIds[1], 1);
    BOOST_CHECK_EQUAL(vCellIds[2], 3);
    BOOST_CHECK_EQUAL(vCellIds[3], 2);
    BOOST_CHECK_EQUAL(vCellIds[4], 4);
    BOOST_CHECK_EQUAL(vCellIds[5], 5);
    BOOST_CHECK_EQUAL(vCellIds[6], 7);
    BOOST_CHECK_EQUAL(vCellIds[7], 6);
}

BOOST_AUTO_TEST_CASE(VoronoiTriangle)
{
    IntegrationTypeTensorProduct<2> quadType1(1, eIntegrationMethod::GAUSS);
    IntegrationTypeTensorProduct<2> quadType2(2, eIntegrationMethod::GAUSS);
    IntegrationTypeTetrahedron tetrahedronType(1);
    BOOST_CHECK_THROW(Visualize::VoronoiGeometryTriangle(quadType1), Exception);
    BOOST_CHECK_THROW(Visualize::VoronoiGeometryTriangle(quadType2), Exception);
    BOOST_CHECK_THROW(Visualize::VoronoiGeometryTriangle(tetrahedronType), Exception);

    IntegrationTypeTriangle triangleTypeFourthOrder(4);
    /*
     *           x
     *          / \
     *         x   x
     *        /\   /\
     *       /   x   \
     *      /    |    \
     *     /     x     \
     *    /    ,´ `,    \
     *   x----x     x----x
     *  /     |     |     \
     * x______x_____x______x
     *
     *  Contains 6 cells and 13 points.
     *
     */
    auto geometry = Visualize::VoronoiGeometryTriangle(triangleTypeFourthOrder);
    std::cout << geometry.voronoiCells.size() << std::endl;
    std::cout << geometry.pointCoordinates.size() << std::endl;
    BOOST_CHECK_EQUAL(geometry.voronoiCells.size(), 6);
    BOOST_CHECK_EQUAL(geometry.pointCoordinates.size(), 13);
}
