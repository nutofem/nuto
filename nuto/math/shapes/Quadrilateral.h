#pragma once

#include "Shape.h"

namespace NuTo
{

class Quadrilateral : public Shape
{
public:
    eShape Enum() const override
    {
        return eShape::Quadrilateral;
    }

    bool IsWithinShape(Eigen::VectorXd xi, double e = 1e-8) const override
    {
        double x = xi[0];
        double y = xi[1];
        return (-1.-e < x) && (x < 1.+e) && (-1.-e < y) && (y < 1.+e);
    }
};

} // namespace NuTo
