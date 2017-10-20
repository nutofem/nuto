#include "BoostUnitTest.h"

#include "base/Group.h"

#include "mechanics/mesh/MeshFem.h"
#include "mechanics/mesh/MeshFemDofConvert.h"

#include "mechanics/interpolation/InterpolationQuadLinear.h"

#include "mechanics/constitutive/laws/LinearElastic.h"

#include "mechanics/integrands/MomentumBalance.h"

#include "mechanics/cell/Cell.h"
#include "mechanics/cell/SimpleAssember.h"

using namespace NuTo;
using namespace NuTo::Groups;


MeshFem QuadPatchTestMesh()
{
    /*
     * Something like this:
     *
     *    3-----------------------2
     * /| | - _        e2       / | -->
     * /| |     -7------------6   | -->
     * /| | e5  /     e4      |   | -->
     * /| |    /              |e1 | -->
     * /| |   /    _____------5   | -->
     * /| |  4-----            \  | --> p
     * /| | /      e0           \ | -->
     * /| |/                     \| -->
     *    0-----------------------1
     *   /_\
     *   ///
     *              (c) ttitsche :)
     */


    MeshFem mesh;
    auto& n0 = mesh.Nodes.Add(Eigen::Vector2d(0, 0));
    auto& n1 = mesh.Nodes.Add(Eigen::Vector2d(10, 0));
    auto& n2 = mesh.Nodes.Add(Eigen::Vector2d(10, 10));
    auto& n3 = mesh.Nodes.Add(Eigen::Vector2d(0, 10));
    
    auto& n4 = mesh.Nodes.Add(Eigen::Vector2d(2, 2));
    auto& n5 = mesh.Nodes.Add(Eigen::Vector2d(8, 3));
    auto& n6 = mesh.Nodes.Add(Eigen::Vector2d(8, 7));
    auto& n7 = mesh.Nodes.Add(Eigen::Vector2d(4, 7));

    const auto& interpolation = mesh.CreateInterpolation(InterpolationQuadLinear(2));

    mesh.Elements.Add({{{n0, n1, n5, n4}, interpolation}});
    mesh.Elements.Add({{{n1, n2, n6, n5}, interpolation}});
    mesh.Elements.Add({{{n7, n6, n2, n3}, interpolation}});
    mesh.Elements.Add({{{n4, n4, n6, n7}, interpolation}});
    mesh.Elements.Add({{{n0, n4, n7, n3}, interpolation}});

    return mesh;
}

BOOST_AUTO_TEST_CASE(PatchTest)
{
    MeshFem mesh = QuadPatchTestMesh(); 
    DofType displ("displacements", 2);
    const auto& interpolation = mesh.CreateInterpolation(InterpolationQuadLinear(2));

    AddDofInterpolation(&mesh, displ, interpolation);
}

