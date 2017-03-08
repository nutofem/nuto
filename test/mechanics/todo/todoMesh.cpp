#include "BoostUnitTest.h"
#include "mechanics/mesh/Mesh.h"

BOOST_AUTO_TEST_CASE(MeshValuePointer)
{
    NuTo::Mesh mesh;
    NuTo::NodeSimple& n0 = mesh.CreateNode(Eigen::Vector2d{1,2});

    for (int i = 0; i < 1e6; ++i)
        mesh.CreateNode(Eigen::Matrix<double, 0, 1>());

    BOOST_CHECK_EQUAL(n0.GetNumValues(), 2);
}
