#include "BoostUnitTest.h"

#include "nuto/mechanics/mesh/MeshFem.h"
#include "nuto/mechanics/mesh/UnitMeshFem.h"
#include "nuto/mechanics/mesh/MeshTopology.h"
#include "nuto/mechanics/mesh/MeshCompanion.h"

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

    MeshTopology meshTopo(&mesh);

    BOOST_CHECK_EQUAL(meshTopo.GetNumEdges(), 20);
    BOOST_CHECK_EQUAL(meshTopo.GetNumFaces(), 11);
    BOOST_CHECK_EQUAL(meshTopo.GetNumVolumes(), 2);
}
