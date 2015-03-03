#ifndef ELEMENTSHAPEFUNCTIONS_H_
#define ELEMENTSHAPEFUNCTIONS_H_


#include <vector>


namespace NuTo
{
namespace ShapeFunctions1D
{
////////////////////////////////////////////////////////////////////////////////////////////////////////////

    void ShapeFunctions1D2N(double rNaturalCoordinates, std::vector<double>& rShapeFunctions);

    void DerivativeShapeFunctions1D2N(double rNaturalCoordinates, std::vector<double>& rDerivativeShapeFunctions);

////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    void ShapeFunctions1D3N(double rNaturalCoordinates, std::vector<double>& rShapeFunctions);

    void DerivativeShapeFunctions1D3N(double rNaturalCoordinates, std::vector<double>& rDerivativeShapeFunctions);

//////////////////////////////////////////////////////////////////////////////////////////////////////////////

    void ShapeFunctions1D4N(double rNaturalCoordinates, std::vector<double>& rShapeFunctions);

    void DerivativeShapeFunctions1D4N(double rNaturalCoordinates, std::vector<double>& rDerivativeShapeFunctions);

//////////////////////////////////////////////////////////////////////////////////////////////////////////////

    void ShapeFunctions1D5N(double rNaturalCoordinates, std::vector<double>& rShapeFunctions);

    void DerivativeShapeFunctions1D5N(double rNaturalCoordinates, std::vector<double>& rDerivativeShapeFunctions);

//////////////////////////////////////////////////////////////////////////////////////////////////////////////

}

namespace ShapeFunctions2D
{
////////////////////////////////////////////////////////////////////////////////////////////////////////////
    void ShapeFunctions2D3N(const double rNaturalCoordinates[2], std::vector<double>& rShapeFunctions);

    void DerivativeShapeFunctions2D3N(const double rNaturalCoordinates[2], std::vector<double>& rDerivativeShapeFunctions);

////////////////////////////////////////////////////////////////////////////////////////////////////////////

    void ShapeFunctions2D4N(const double rNaturalCoordinates[2], std::vector<double>& rShapeFunctions);

    void DerivativeShapeFunctions2D4N(const double rNaturalCoordinates[2], std::vector<double>& rDerivativeShapeFunctions);

////////////////////////////////////////////////////////////////////////////////////////////////////////////

    void ShapeFunctions2D6N(const double rNaturalCoordinates[2], std::vector<double>& rShapeFunctions);

    void DerivativeShapeFunctions2D6N(const double rNaturalCoordinates[2], std::vector<double>& rDerivativeShapeFunctions);

////////////////////////////////////////////////////////////////////////////////////////////////////////////

    void ShapeFunctions2D10N(const double rNaturalCoordinates[2], std::vector<double>& rShapeFunctions);

    void DerivativeShapeFunctions2D10N(const double rNaturalCoordinates[2], std::vector<double>& rDerivativeShapeFunctions);

////////////////////////////////////////////////////////////////////////////////////////////////////////////

    void ShapeFunctions2D15N(const double rNaturalCoordinates[2], std::vector<double>& rShapeFunctions);

    void DerivativeShapeFunctions2D15N(const double rNaturalCoordinates[2], std::vector<double>& rDerivativeShapeFunctions);

////////////////////////////////////////////////////////////////////////////////////////////////////////////

}


namespace ShapeFunctions3D
{

}

}



#endif /* ELEMENTSHAPEFUNCTIONS_H_ */
