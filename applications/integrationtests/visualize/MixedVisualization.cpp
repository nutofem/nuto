#include "BoostUnitTest.h"
#include <eigen3/Eigen/Core>
#include "visualize/VisualizeUnstructuredGrid.h"

BOOST_AUTO_TEST_CASE(TriangleAndQuad)
{
    NuTo::VisualizeUnstructuredGrid visualizer;
  
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
    Eigen::Vector3d p0(0,0,0);
    Eigen::Vector3d p1(1,0,0);
    Eigen::Vector3d p2(1,2,0);
    Eigen::Vector3d p3(0,2,0);
    Eigen::Vector3d p4(2,1,0);

    BOOST_CHECK_EQUAL(visualizer.AddPoint(p0.data()), 0);
    BOOST_CHECK_EQUAL(visualizer.AddPoint(p1.data()), 1);
    BOOST_CHECK_EQUAL(visualizer.AddPoint(p2.data()), 2);
    BOOST_CHECK_EQUAL(visualizer.AddPoint(p3.data()), 3);
    BOOST_CHECK_EQUAL(visualizer.AddPoint(p4.data()), 4);

    std::vector<unsigned int> quadIds = {0,1,2,3};
    std::vector<unsigned int> triangleIds = {1,4,2};

    visualizer.AddQuadCell(quadIds.data());
    visualizer.AddTriangleCell(triangleIds.data());

    BOOST_CHECK_NO_THROW(visualizer.ExportVtuDataFile("MixedVisualization.vtu"));
}
