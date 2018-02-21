#include "BoostUnitTest.h"
#include "visualize/Visualizer.h"
#include "visualize/PointHandler.h"
#include "mechanics/integrationtypes/IntegrationTypeTensorProduct.h"

#include "TestStructure.h"

using namespace NuTo;
using namespace NuTo::Visualize;

BOOST_AUTO_TEST_CASE(PointVisualizer)
{
    NuTo::Test::VisualizeTestStructure s;
    auto cells = s.Cells();

    IntegrationTypeTensorProduct<2> integrationType(2, eIntegrationMethod::GAUSS);

    std::string filename = "PointHandlerOutput.vtu";

    Visualizer visualize(cells, PointHandler(integrationType));
    visualize.WriteVtuFile(filename, false);
    UnstructuredGridCheck::CheckNum(filename, 8, 0);
}
