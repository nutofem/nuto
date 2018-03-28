#include "nuto/visualize/AverageHandler.h"
#include "nuto/visualize/Visualizer.h"
#include "nuto/visualize/XMLWriter.h"
#include "TestStructure.h"

using namespace NuTo;
using namespace NuTo::Visualize;

BOOST_AUTO_TEST_CASE(GroupAverage)
{
    NuTo::Test::VisualizeTestStructure s;
    auto cells = s.Cells();

    std::string filename = "AverageOutput.vtu";
    Visualizer visualize(cells, AverageHandler());
    visualize.WriteVtuFile(filename, false);
    UnstructuredGridCheck::CheckNum(filename, 8, 2);
}
