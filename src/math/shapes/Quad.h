#pragma once

#include "Shape.h"

namespace NuTo
{

class Quad : public Shape
{
public:
    eShape Enum() const override
    {
        return eShape::Quadrilateral;
    }

    bool IsWithinShape(Eigen::VectorXd xi) const override
    {
        double x = xi[0];
        double y = xi[1];
        return (-1.<x)&&(x<1.)&&(-1.<y)&&(y<1.);
    }
};

} // namespace NuTo
