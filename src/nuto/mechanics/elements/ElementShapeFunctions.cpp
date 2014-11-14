
#include "nuto/mechanics/elements/ElementShapeFunctions.h"
#include <assert.h>


namespace NuTo
{
namespace ShapeFunctions1D
{
////////////////////////////////////////////////////////////////////////////////////////////////////////////
    void ShapeFunctions1D2N(double rLocalCoordinates, std::vector<double>& rShapeFunctions)
    {
        assert(rShapeFunctions.size()==2);
        rShapeFunctions[0] = 0.5*(1.-rLocalCoordinates);
        rShapeFunctions[1] = 0.5*(1.+rLocalCoordinates);
    }

    void DerivativeShapeFunctions1D2N(double rLocalCoordinates, std::vector<double>& rDerivativeShapeFunctions)
    {
        assert(rDerivativeShapeFunctions.size()==2);
        rDerivativeShapeFunctions[0] = -0.5;
        rDerivativeShapeFunctions[1] = 0.5;
    }

////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    void ShapeFunctions1D3N(double rLocalCoordinates, std::vector<double>& rShapeFunctions)
    {
        assert(rShapeFunctions.size()==3);
        rShapeFunctions[0] = 0.5*(1.-rLocalCoordinates)-0.5*(1.-rLocalCoordinates*rLocalCoordinates);
        rShapeFunctions[1] = 1.-rLocalCoordinates*rLocalCoordinates;
        rShapeFunctions[2] = 0.5*(1.+rLocalCoordinates)-0.5*(1.-rLocalCoordinates*rLocalCoordinates);
    }

    void DerivativeShapeFunctions1D3N(double rLocalCoordinates, std::vector<double>& rDerivativeShapeFunctions)
    {
        assert(rDerivativeShapeFunctions.size()==3);
        rDerivativeShapeFunctions[0] = -0.5 + rLocalCoordinates;
        rDerivativeShapeFunctions[1] = -2.0 * rLocalCoordinates;
        rDerivativeShapeFunctions[2] =  0.5 + rLocalCoordinates;
    }

//////////////////////////////////////////////////////////////////////////////////////////////////////////////


}

namespace ShapeFunctions2D
{

}


namespace ShapeFunctions3D
{

}

}
