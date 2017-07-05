#include "BoostUnitTest.h"
#include <eigen3/Eigen/Core>
#include "visualize/UnstructuredGrid.h"

BOOST_AUTO_TEST_CASE(TriangleAndQuad)
{
    NuTo::Visualize::UnstructuredGrid visualizer;

    /**
     * 3_____2
     * |     |\
     * |     | \4
     * |     | /
     * 0_____1/
     *
     *
     *
     */
    Eigen::Vector3d p0(0, 0, 0);
    Eigen::Vector3d p1(1, 0, 0);
    Eigen::Vector3d p2(1, 2, 0);
    Eigen::Vector3d p3(0, 2, 0);
    Eigen::Vector3d p4(2, 1, 0);

    BOOST_CHECK_EQUAL(visualizer.AddPoint(p0), 0);
    BOOST_CHECK_EQUAL(visualizer.AddPoint(p1), 1);
    BOOST_CHECK_EQUAL(visualizer.AddPoint(p2), 2);
    BOOST_CHECK_EQUAL(visualizer.AddPoint(p3), 3);
    BOOST_CHECK_EQUAL(visualizer.AddPoint(p4), 4);

    std::vector<int> quadIds = {0, 1, 2, 3};
    std::vector<int> triangleIds = {1, 4, 2};

    visualizer.AddCell(quadIds, NuTo::eCellTypes::QUAD);
    visualizer.AddCell(triangleIds, NuTo::eCellTypes::TRIANGLE);

    BOOST_CHECK_NO_THROW(visualizer.ExportVtuDataFile("MixedVisualization.vtu"));
}
