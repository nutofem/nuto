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

    bool IsWithinShape(Eigen::VectorXd xi, double e = 1e-8) const override
    {
        double x = xi[0];
        double y = xi[1];
        return (0.-e < x) && (x < 1.+e) && (0.-e < y) && (y < 1.+e) && (x+y < 1.+e);
    }
};

} // namespace NuTo
