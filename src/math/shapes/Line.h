#pragma once

#include "Shape.h"

namespace NuTo
{

class Line : public Shape
{
public:
    eShape Enum() const override
    {
        return eShape::Line;
    }

    bool IsWithinShape(Eigen::VectorXd xi) const override
    {
        double x = xi[0];
        return (-1.<x)&&(x<1.);
    }
};

} // namespace NuTo
