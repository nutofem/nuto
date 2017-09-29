#pragma once
#include <eigen3/Eigen/Dense>
#include "mechanics/interpolation/TypeDefs.h"

namespace NuTo
{

class ElementInterface
{
public:
    virtual ~ElementInterface() = default;

    //! @brief extracts all node values of this element
    virtual NodeValues ExtractNodeValues() const = 0;
    virtual int GetDofDimension() const = 0;
    virtual int GetNumNodes() const = 0;
    virtual NMatrix GetNMatrix(NaturalCoords ipCoords) const = 0;
    virtual ShapeFunctions GetShapeFunctions(NaturalCoords ipCoords) const = 0;
    virtual DerivativeShapeFunctionsNatural GetDerivativeShapeFunctions(NaturalCoords ipCoords) const = 0;
};

Eigen::VectorXd Interpolate(const ElementInterface& element, NaturalCoords ipCoords)
{
    return element.GetNMatrix(ipCoords) * element.ExtractNodeValues();
}

} /* NuTo */
