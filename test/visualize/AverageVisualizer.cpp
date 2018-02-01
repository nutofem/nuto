#include "BoostUnitTest.h"
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/xml_parser.hpp>
#include "visualize/AverageHandler.h"
#include "visualize/AverageGeometries.h"
#include "visualize/Visualizer.h"
#include "visualize/XMLWriter.h"
#include "TestStructure.h"

namespace pt = boost::property_tree;
using namespace NuTo;

struct UnstructuredGridCheck
{
public:
    static void Average(std::string filename)
    {
        pt::ptree tree;
        pt::read_xml(filename, tree);
        int numOfPoints = tree.get<int>("VTKFile.UnstructuredGrid.Piece.<xmlattr>.NumberOfPoints");
        int numOfCells = tree.get<int>("VTKFile.UnstructuredGrid.Piece.<xmlattr>.NumberOfCells");
        BOOST_CHECK_EQUAL(numOfPoints, 8);
        BOOST_CHECK_EQUAL(numOfCells, 2);
    }
};


BOOST_AUTO_TEST_CASE(GroupAverage)
{
    using namespace Visualize;
    NuTo::DofType dof("NodeCoordinatesDiv10", 2);
    NuTo::Test::VisualizeTestStructure s(dof);
    auto cells = s.Cells();

    std::string filename = "AverageOutput.vtu";
    Visualizer<AverageHandler> visualize(cells, AverageGeometryQuad());
    visualize.WriteVtuFile(filename, false);
    UnstructuredGridCheck::Average(filename);
}
