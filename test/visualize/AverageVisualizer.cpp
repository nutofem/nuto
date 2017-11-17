#include "BoostUnitTest.h"
#include "visualize/AverageVisualizer.h"
#include "visualize/XMLWriter.h"
#include "TestStructure.h"

struct UnstructuredGridCheck
{
public:
    static void Average(const NuTo::Visualize::UnstructuredGrid& grid)
    {
        BOOST_CHECK_EQUAL(grid.mPoints.size(), 8);
        BOOST_CHECK_EQUAL(grid.mCells.size(), 2);
    }
};


BOOST_AUTO_TEST_CASE(GroupAverage)
{
    NuTo::DofType dof("NodeCoordinatesDiv10", 2);
    NuTo::Test::VisualizeTestStructure s(dof);


    NuTo::Visualize::CellGeometry cellGeometry;
    NuTo::InterpolationQuadLinear interpolation(2);
    cellGeometry.mCornerCoords = {interpolation.GetLocalCoords(0), interpolation.GetLocalCoords(1),
                                  interpolation.GetLocalCoords(2), interpolation.GetLocalCoords(3)};
    cellGeometry.mCellType = NuTo::eCellTypes::QUAD;


    NuTo::Visualize::GroupAverage v(s.Cells(), cellGeometry);
    v.GeometryToGrid();
    v.DofsToGrid({dof});
    v.IpValuesToGrid();

    const auto& grid = v.GetUnstructuredGrid();

    UnstructuredGridCheck::Average(grid);

    NuTo::Visualize::XMLWriter::Export("PdeAverage.vtu", grid, false);
}
