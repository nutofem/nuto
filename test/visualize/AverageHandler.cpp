#include "visualize/AverageGeometries.h"
#include "visualize/Visualizer.h"
#include "visualize/XMLWriter.h"
#include "TestStructure.h"

using namespace NuTo;
using namespace NuTo::Visualize;

BOOST_AUTO_TEST_CASE(GroupAverage)
{
    NuTo::Test::VisualizeTestStructure s;
    auto cells = s.Cells();

    std::string filename = "AverageOutput.vtu";
    Visualizer visualize(cells, AverageHandler(AverageGeometryQuad()));
    visualize.WriteVtuFile(filename, false);
    UnstructuredGridCheck::CheckNum(filename, 8, 2);
}
