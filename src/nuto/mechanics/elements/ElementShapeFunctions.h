#ifndef ELEMENTSHAPEFUNCTIONS_H_
#define ELEMENTSHAPEFUNCTIONS_H_


#include <vector>


namespace NuTo
{
namespace ShapeFunctions1D
{
////////////////////////////////////////////////////////////////////////////////////////////////////////////
    void ShapeFunctions1D2N(double rLocalCoordinates, std::vector<double>& rShapeFunctions);

    void DerivativeShapeFunctions1D2N(double rLocalCoordinates, std::vector<double>& rDerivativeShapeFunctions);

////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    void ShapeFunctions1D3N(double rLocalCoordinates, std::vector<double>& rShapeFunctions);

    void DerivativeShapeFunctions1D3N(double rLocalCoordinates, std::vector<double>& rDerivativeShapeFunctions);

//////////////////////////////////////////////////////////////////////////////////////////////////////////////


}

namespace ShapeFunctions2D
{

}


namespace ShapeFunctions3D
{

}

}



#endif /* ELEMENTSHAPEFUNCTIONS_H_ */
