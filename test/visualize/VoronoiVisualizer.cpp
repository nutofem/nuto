#include "BoostUnitTest.h"
#include "BoostUnitTest.h"
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/xml_parser.hpp>
#include "visualize/TensorProductVoronoiHandler.h"
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
    NuTo::DofType dof("NodeCoordinatesDiv10", 2);
    NuTo::Test::VisualizeTestStructure s(dof);
    auto cells = s.Cells();

    std::string filename = "VoronoiOutput.vtu";
    Visualize::Visualizer<Visualize::TensorProductVoronoiHandler<2>> visualize(cells, 2);
    visualize.WriteVtuFile(filename);
    UnstructuredGridCheck::Voronoi(filename);
}
