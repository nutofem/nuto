#pragma once
#include <eigen3/Eigen/Core>
#include "nuto/mechanics/interpolation/TypeDefs.h"

namespace NuTo
{

class ElementInterface
{
public:
    virtual ~ElementInterface() noexcept = default;

    virtual int GetNumNodes() const = 0;
    virtual int GetDofDimension() const = 0;
    virtual Eigen::MatrixXd GetNMatrix(NaturalCoords ipCoords) const = 0;
    virtual Eigen::VectorXd GetShapeFunctions(NaturalCoords ipCoords) const = 0;
    virtual Eigen::MatrixXd GetDerivativeShapeFunctions(NaturalCoords ipCoords) const = 0;
};

} /* NuTo */
