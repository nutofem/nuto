#pragma once
#include<eigen3/Eigen/Dense>
#include "TypeDefs.h"

namespace NuTo
{

class ElementInterpolationBase
{
public:
    ElementInterpolationBase(){}

    //! @brief extracts all node values of this element
    //! @remark virtual to make it testable
    virtual NodeValues ExtractNodeValues() const = 0;

    virtual int GetDofDimension() const = 0;

    virtual int GetNumNodes() const = 0;

    virtual Eigen::VectorXd GetShapeFunctions(Eigen::VectorXd ipCoords) const = 0;

    virtual Eigen::MatrixXd GetDerivativeShapeFunctions(Eigen::VectorXd ipCoords) const = 0;
};
} /* NuTo */