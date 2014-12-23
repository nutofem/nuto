
#include "nuto/mechanics/elements/ElementShapeFunctions.h"
#include <assert.h>


namespace NuTo
{
namespace ShapeFunctions1D // interval -1 to 1
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

    void ShapeFunctions1D4N(double rLocalCoordinates, std::vector<double>& rShapeFunctions)
    {
        assert(rShapeFunctions.size()==4);
        double x = rLocalCoordinates;
        rShapeFunctions[0] = -  9./16. * (x+1./3) * (x-1./3) * (x-1);
        rShapeFunctions[1] =   27./16. * (x+1.)   * (x-1./3) * (x-1);
        rShapeFunctions[2] = - 27./16  * (x+1.)   * (x+1./3) * (x-1);
        rShapeFunctions[3] =    9./16. * (x-1./3) * (x+1./3) * (x+1);
    }

    void DerivativeShapeFunctions1D4N(double rLocalCoordinates, std::vector<double>& rDerivativeShapeFunctions)
    {
        assert(rDerivativeShapeFunctions.size()==4);
        double x = rLocalCoordinates;
        rDerivativeShapeFunctions[0] = (-27.0*x*x + 18.0*x +  1.0) / 16.;
        rDerivativeShapeFunctions[1] = ( 81.0*x*x - 18.0*x - 27.0) / 16.;
        rDerivativeShapeFunctions[2] = (-81.0*x*x - 18.0*x + 27.0) / 16.;
        rDerivativeShapeFunctions[3] = ( 27.0*x*x + 18.0*x -  1.0) / 16.;
    }
}

namespace ShapeFunctions2D
{

}


namespace ShapeFunctions3D
{

}

}
