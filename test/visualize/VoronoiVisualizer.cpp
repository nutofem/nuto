#include "BoostUnitTest.h"
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/xml_parser.hpp>
#include "visualize/VoronoiHandler.h"
#include "visualize/VoronoiGeometries.h"
#include "visualize/Visualizer.h"
#include "visualize/XMLWriter.h"
#include "TestStructure.h"

namespace pt = boost::property_tree;
using namespace NuTo;

struct UnstructuredGridCheck
{
public:
    static void Voronoi(std::string filename)
    {
        pt::ptree tree;
        pt::read_xml(filename, tree);
        int numOfPoints = tree.get<int>("VTKFile.UnstructuredGrid.Piece.<xmlattr>.NumberOfPoints");
        int numOfCells = tree.get<int>("VTKFile.UnstructuredGrid.Piece.<xmlattr>.NumberOfCells");
        BOOST_CHECK_EQUAL(numOfPoints, 18);
        BOOST_CHECK_EQUAL(numOfCells, 8);
    }
};

BOOST_AUTO_TEST_CASE(GroupVoronoiTensorProduct2D)
{
    using namespace NuTo::Visualize;
    NuTo::DofType dof("NodeCoordinatesDiv10", 2);
    NuTo::Test::VisualizeTestStructure s(dof);
    auto cells = s.Cells();

    std::string filename = "VoronoiOutput.vtu";
    Visualizer<VoronoiHandler> visualize(cells, VoronoiGeometryQuad(2));
    visualize.WriteVtuFile(filename, false);
    UnstructuredGridCheck::Voronoi(filename);
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
