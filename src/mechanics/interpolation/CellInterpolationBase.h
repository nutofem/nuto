#pragma once
#include <eigen3/Eigen/Dense>
#include "TypeDefs.h"

namespace NuTo
{

class CellInterpolationBase
{
public:
    virtual ~CellInterpolationBase() = default;

    //! @brief extracts all node values of this element
    virtual NodeValues ExtractNodeValues() const = 0;
    virtual int GetDofDimension() const = 0;
    virtual int GetNumNodes() const = 0;
    virtual NMatrix GetNMatrix(NaturalCoords ipCoords) const = 0;
    virtual ShapeFunctions GetShapeFunctions(NaturalCoords ipCoords) const = 0;
    virtual DerivativeShapeFunctionsNatural GetDerivativeShapeFunctions(NaturalCoords ipCoords) const = 0;
};
} /* NuTo */
