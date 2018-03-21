#include "BoostUnitTest.h"
#include "nuto/visualize/AverageHandler.h"
#include "nuto/visualize/Visualizer.h"
#include "TestStructure.h"

using namespace NuTo;
using namespace NuTo::Visualize;


BOOST_AUTO_TEST_CASE(EmptyCells)
{
    Visualizer visualize({}, AverageHandler());

    NuTo::Test::VisualizeTestStructure s;
    auto cells = s.Cells();

    visualize = Visualizer(cells, AverageHandler());
    std::string filename = "EmptyCellsOutput.vtu";
    visualize.WriteVtuFile(filename, false);
    UnstructuredGridCheck::CheckNum(filename, 8, 2);
}

BOOST_AUTO_TEST_CASE(ReplacingCells)
{
    NuTo::Test::VisualizeTestStructure s;
    auto cells = s.Cells();

    Visualizer visualize(cells, AverageHandler());
    visualize = Visualizer(cells, AverageHandler());
    std::string filename = "ReplacingCellsOutput.vtu";
    visualize.WriteVtuFile(filename, false);
    UnstructuredGridCheck::CheckNum(filename, 8, 2);
}
