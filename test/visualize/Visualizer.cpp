#include "BoostUnitTest.h"
#include "visualize/AverageGeometries.h"
#include "visualize/Visualizer.h"
#include "TestStructure.h"

using namespace NuTo;
using namespace NuTo::Visualize;

Visualizer CreateVisualizer(Group<CellInterface> cells)
{
    AverageHandler handler(AverageGeometryQuad());
    Visualizer visu(cells, handler);
    return visu; // tests the _copy/clone_, because `handler` is out of scope now.
}

BOOST_AUTO_TEST_CASE(HandlerPolymorphism)
{
    NuTo::DofType dof("NodeCoordinatesDiv10", 2);
    NuTo::Test::VisualizeTestStructure s(dof);
    auto cells = s.Cells();

    std::string filename = "PolyOutput.vtu";
    Visualizer visualize = CreateVisualizer(cells);
    visualize.WriteVtuFile(filename, false);
    UnstructuredGridCheck::CheckNum(filename, 8, 2);
}


BOOST_AUTO_TEST_CASE(EmptyCells)
{
    Visualizer visualize({}, AverageGeometryQuad());

    NuTo::DofType dof("NodeCoordinatesDiv10", 2);
    NuTo::Test::VisualizeTestStructure s(dof);
    auto cells = s.Cells();

    visualize.SetCells(cells);
    std::string filename = "EmptyCellsOutput.vtu";
    visualize.WriteVtuFile(filename, false);
    UnstructuredGridCheck::CheckNum(filename, 8, 2);
}

BOOST_AUTO_TEST_CASE(ReplacingCells)
{

    NuTo::DofType dof("NodeCoordinatesDiv10", 2);
    NuTo::Test::VisualizeTestStructure s(dof);
    auto cells = s.Cells();

    Visualizer visualize = CreateVisualizer(cells);
    visualize.SetCells(cells);
    std::string filename = "ReplacingCellsOutput.vtu";
    visualize.WriteVtuFile(filename, false);
    UnstructuredGridCheck::CheckNum(filename, 8, 2);
}
