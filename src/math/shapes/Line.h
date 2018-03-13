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

    bool IsWithinShape(Eigen::VectorXd xi, double e = 1e-8) const override
    {
        double x = xi[0];
        return (-1.-e < x) && (x < 1.+e);
    }
};

} // namespace NuTo
