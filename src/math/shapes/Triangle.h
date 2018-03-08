#pragma once

#include "Shape.h"

namespace NuTo
{

class Triangle : public Shape
{
public:
    eShape Enum() const override
    {
        return eShape::Triangle;
    }

    bool IsWithinShape(Eigen::VectorXd xi) const override
    {
        double x = xi[0];
        double y = xi[1];
        return (0.<x)&&(x<1.)&&(0.<y)&&(y<1.)&&(x+y<1.);
    }
};

} // namespace NuTo
