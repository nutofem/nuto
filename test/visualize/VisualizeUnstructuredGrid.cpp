#include "BoostUnitTest.h"
#include "visualize/VisualizeUnstructuredGrid.h"
#include "visualize/VisualizeException.h"

BOOST_AUTO_TEST_CASE(DefinitionOrder)
{
    NuTo::VisualizeUnstructuredGrid visu;
    visu.DefinePointData("Stuff");
    visu.DefineCellData("Stuff");
    int pointId = visu.AddPoint(Eigen::Vector3d::Zero());
    NuTo::CellBase cell({pointId}, 1, NuTo::eCellTypes::VERTEX);
    int cellId = visu.AddCell(cell);

    // add data to undefined data field
    BOOST_CHECK_THROW(visu.SetPointData(pointId, "OtherStuff", 42.), NuTo::VisualizeException);
    BOOST_CHECK_THROW(visu.SetCellData(cellId, "OtherStuff", 42.), NuTo::VisualizeException);

    // definine data after adding points/cells
    BOOST_CHECK_THROW(visu.DefinePointData("OtherStuff"), NuTo::VisualizeException);
    BOOST_CHECK_THROW(visu.DefineCellData("OtherStuff"), NuTo::VisualizeException);
}
