#include "nuto/visualize/VoronoiGeometries.h"
#include "nuto/visualize/Visualizer.h"
#include "nuto/visualize/XMLWriter.h"
#include "TestStructure.h"

using namespace NuTo;


BOOST_AUTO_TEST_CASE(GroupVoronoiTensorProduct2D)
{
    using namespace NuTo::Visualize;
    NuTo::Test::VisualizeTestStructure s;
    auto cells = s.Cells();

    std::string filename = "VoronoiOutput.vtu";
    Visualizer visualize(cells, VoronoiHandler(VoronoiGeometryQuad(2)));
    visualize.WriteVtuFile(filename, false);
    UnstructuredGridCheck::CheckNum(filename, 18, 8);
}
