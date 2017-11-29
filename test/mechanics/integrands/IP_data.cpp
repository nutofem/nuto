#include "BoostUnitTest.h"


#include "mechanics/nodes/NodeSimple.h"
#include "mechanics/interpolation/InterpolationTrussLinear.h"
#include "mechanics/mesh/MeshFem.h"

using namespace NuTo;


MeshFem Mesh1D()
{
    MeshFem mesh;
    const InterpolationSimple& interpolation = mesh.CreateInterpolation(InterpolationTrussLinear(1));
    NodeSimple* nr = nullptr;
    for (unsigned int i = 0; i < 21; ++i)
    {
        NodeSimple& nl = mesh.Nodes.Add({i * 0.5});
        if (i > 0)
            mesh.Elements.Add({{{nl, *nr}, interpolation}});
        nr = &nl;
    }


    return mesh;
}

BOOST_AUTO_TEST_CASE(IP_data)
{
    MeshFem mesh = Mesh1D();
    DofType displ("displacements", 1);
    const auto& interpolation = mesh.CreateInterpolation(InterpolationTrussLinear(1));
}
