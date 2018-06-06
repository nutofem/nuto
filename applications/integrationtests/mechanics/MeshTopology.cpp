#include "BoostUnitTest.h"

#include "nuto/mechanics/mesh/MeshFem.h"
#include "nuto/mechanics/mesh/UnitMeshFem.h"
#include "nuto/mechanics/mesh/MeshTopology.h"
#include "nuto/mechanics/mesh/MeshCompanion.h"
#include "nuto/mechanics/interpolation/InterpolationBrickLinear.h"
#include "nuto/mechanics/interpolation/InterpolationQuadLinear.h"
#include "nuto/mechanics/interpolation/InterpolationTrussLinear.h"

using namespace NuTo;

/* Test mesh consisting of two Hexes with
 * all edges and faces (all linear)
 *
 *      x-----x-----x
 *     /|    /     /|
 *    x-|---x-----x |
 *    | |   |     | |
 *    | x---|-x---|-x
 *    |/    |     |/
 *    x-----x-----x
*/
BOOST_AUTO_TEST_CASE(CreateMeshTopology)
{
    MeshFem mesh = UnitMeshFem::CreateBricks(2, 1, 1);
    auto allVolumeElements = mesh.ElementsTotal();
    AddEdgeElements(&mesh, allVolumeElements);
    AddFaceElements(&mesh, allVolumeElements);

    MeshTopology meshTopo(mesh.ElementsTotal());

    BOOST_CHECK_EQUAL(meshTopo.GetNumEdges(), 20);
    BOOST_CHECK_EQUAL(meshTopo.GetNumFaces(), 11);
    BOOST_CHECK_EQUAL(meshTopo.GetNumVolumes(), 2);
}

/* Test mesh consisting of two Hexes with
 *
 *    Nodes:              Edges:          Faces, Volumes
 *
 *      7-----6-----11                         6
 *     /|    /     /|                 face0   /|
 *    3-|---2-----9 |            2           2 |
 *    | |   |     | |            |      vol0 | |  vol1
 *    | 4---|-5---|-10           e1 5        | 5
 *    |/    |     |/             | / e2      |/
 *    0-----1-----8        0-e0->1           1
*/
BOOST_AUTO_TEST_CASE(Adjacency)
{
    NodeSimple nd0(Eigen::Vector3d(0., 0., 0.));
    NodeSimple nd1(Eigen::Vector3d(1., 0., 0.));
    NodeSimple nd2(Eigen::Vector3d(1., 1., 0.));
    NodeSimple nd3(Eigen::Vector3d(0., 1., 0.));
    NodeSimple nd4(Eigen::Vector3d(0., 0., 1.));
    NodeSimple nd5(Eigen::Vector3d(1., 0., 1.));
    NodeSimple nd6(Eigen::Vector3d(1., 1., 1.));
    NodeSimple nd7(Eigen::Vector3d(0., 1., 1.));

    NodeSimple nd8(Eigen::Vector3d(2., 0., 0.));
    NodeSimple nd9(Eigen::Vector3d(2., 1., 0.));
    NodeSimple nd10(Eigen::Vector3d(2., 0., 1.));
    NodeSimple nd11(Eigen::Vector3d(2., 1., 1.));

    InterpolationBrickLinear ipol3d;
    InterpolationQuadLinear ipol2d;
    InterpolationTrussLinear ipol1d;

    MeshFem mesh;
    auto& vol0 = mesh.Elements.Add({{{nd0, nd1, nd2, nd3, nd4, nd5, nd6, nd7}, ipol3d}});
    auto& vol1 = mesh.Elements.Add({{{nd1, nd8, nd9, nd2, nd5, nd10, nd11, nd6}, ipol3d}});

    auto& edge0 = mesh.Elements.Add({{{nd0, nd1}, ipol1d}});
    auto& edge1 = mesh.Elements.Add({{{nd1, nd2}, ipol1d}});
    auto& edge2 = mesh.Elements.Add({{{nd1, nd5}, ipol1d}});

    auto& face0 = mesh.Elements.Add({{{nd1, nd5, nd6, nd2}, ipol2d}});

    MeshTopology meshTopo(mesh.ElementsTotal());

    auto edgesOfFace0 = meshTopo.GetAdjacentEdges(face0);
    BOOST_CHECK(edgesOfFace0.Contains(edge1));
    BOOST_CHECK_EQUAL(edgesOfFace0.Size(), 2);

    auto volumesOfFace0 = meshTopo.GetAdjacentVolumes(face0);
    BOOST_CHECK_EQUAL(volumesOfFace0.Size(), 2);
    BOOST_CHECK(volumesOfFace0.Contains(vol0));
    BOOST_CHECK(volumesOfFace0.Contains(vol1));

    auto facesOfFace0 = meshTopo.GetAdjacentFaces(face0);
    BOOST_CHECK_EQUAL(facesOfFace0.Size(), 1);
    BOOST_CHECK(facesOfFace0.Contains(face0));
}
