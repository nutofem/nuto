#include "BoostUnitTest.h"
#include "nuto/visualize/Visualizer.h"
#include "nuto/visualize/PointHandler.h"
#include "nuto/mechanics/integrationtypes/IntegrationTypeTensorProduct.h"

#include "TestStructure.h"

using namespace NuTo;
using namespace NuTo::Visualize;

BOOST_AUTO_TEST_CASE(PointVisualizer)
{
    NuTo::Test::VisualizeTestStructure s;
    auto cells = s.Cells();

    Eigen::Vector2d a(0.5, 0.5);
    std::string filename = "PointHandlerOutput.vtu";

    Visualizer visualize(cells, PointHandler({a, a, a, a}));
    visualize.WriteVtuFile(filename, false);
    UnstructuredGridCheck::CheckNum(filename, 8, 0);
}
