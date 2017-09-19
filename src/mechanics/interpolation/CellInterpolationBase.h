#pragma once
#include <eigen3/Eigen/Dense>
#include "TypeDefs.h"

namespace NuTo
{

class CellInterpolationBase
{
public:
    CellInterpolationBase()
    {
    }

    //! @brief extracts all node values of this element
    virtual NodeValues ExtractNodeValues() const = 0;

    virtual int GetDofDimension() const = 0;

    virtual int GetNumNodes() const = 0;

    virtual Eigen::VectorXd GetShapeFunctions(Eigen::VectorXd ipCoords) const = 0;

    virtual Eigen::MatrixXd GetDerivativeShapeFunctions(Eigen::VectorXd ipCoords) const = 0;
};
} /* NuTo */
