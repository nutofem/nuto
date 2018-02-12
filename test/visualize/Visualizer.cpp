#include "BoostUnitTest.h"
#include "visualize/AverageGeometries.h"
#include "visualize/Visualizer.h"
#include "TestStructure.h"

using namespace NuTo;
using namespace NuTo::Visualize;


BOOST_AUTO_TEST_CASE(EmptyCells)
{
    Visualizer visualize({}, AverageHandler(AverageGeometryQuad()));

    NuTo::DofType dof("NodeCoordinatesDiv10", 2);
    NuTo::Test::VisualizeTestStructure s(dof);
    auto cells = s.Cells();

    visualize = Visualizer(cells, AverageHandler(AverageGeometryQuad()));
    std::string filename = "EmptyCellsOutput.vtu";
    visualize.WriteVtuFile(filename, false);
    UnstructuredGridCheck::CheckNum(filename, 8, 2);
}

BOOST_AUTO_TEST_CASE(ReplacingCells)
{

    NuTo::DofType dof("NodeCoordinatesDiv10", 2);
    NuTo::Test::VisualizeTestStructure s(dof);
    auto cells = s.Cells();

    Visualizer visualize(cells, AverageHandler(AverageGeometryQuad()));
    visualize = Visualizer(cells, AverageHandler(AverageGeometryQuad()));
    std::string filename = "ReplacingCellsOutput.vtu";
    visualize.WriteVtuFile(filename, false);
    UnstructuredGridCheck::CheckNum(filename, 8, 2);
}
