#include "BoostUnitTest.h"
#include "visualize/VoronoiVisualizer.h"
#include "visualize/XMLWriter.h"
#include "TestStructure.h"
#include "visualize/CellGeometryVoronoi.h"
#include "TensorProductGeometryVoronoi.h"

struct UnstructuredGridCheck
{
public:
    static void Voronoi(const NuTo::Visualize::UnstructuredGrid& grid)
    {
        BOOST_CHECK_EQUAL(grid.mPoints.size(), 16 * 2);
        BOOST_CHECK_EQUAL(grid.mCells.size(), 8);

        BoostUnitTest::CheckEigenMatrix(grid.mCells[0].GetData(0), Eigen::Vector3d(1, 1, 1));
        BoostUnitTest::CheckEigenMatrix(grid.mCells[0].GetData(1), Eigen::Vector3d(10, 10, 10));

        BoostUnitTest::CheckEigenMatrix(grid.mCells[1].GetData(0), Eigen::Vector3d(2, 2, 2));
        BoostUnitTest::CheckEigenMatrix(grid.mCells[1].GetData(1), Eigen::Vector3d(20, 20, 20));

        BoostUnitTest::CheckEigenMatrix(grid.mCells[2].GetData(0), Eigen::Vector3d(3, 3, 3));
        BoostUnitTest::CheckEigenMatrix(grid.mCells[2].GetData(1), Eigen::Vector3d(30, 30, 30));

        BoostUnitTest::CheckEigenMatrix(grid.mCells[3].GetData(0), Eigen::Vector3d(4, 4, 4));
        BoostUnitTest::CheckEigenMatrix(grid.mCells[3].GetData(1), Eigen::Vector3d(40, 40, 40));
    }
};

BOOST_AUTO_TEST_CASE(GroupVoronoiTensorProduct2D)
{
    NuTo::DofType dof("NodeCoordinatesDiv10", 2);
    NuTo::Test::VisualizeTestStructure s(dof);

    NuTo::Visualize::CellGeometryVoronoi cellGeometry;

    std::vector<double> points = {-0.5, 0.5};
    NuTo::Test::TensorProductGeometryVoronoi<2>::buildCellGeometryVoronoi(points, cellGeometry);

    NuTo::Visualize::GroupVoronoi v(s.Cells(), cellGeometry);
    v.GeometryToGrid();
    v.DofsToGrid({dof});
    v.IpValuesToGrid();

    const auto& grid = v.GetUnstructuredGrid();

    UnstructuredGridCheck::Voronoi(grid);

    NuTo::Visualize::XMLWriter::Export("PdeTensorVoronoiTest2D.vtu", grid, false);
};
