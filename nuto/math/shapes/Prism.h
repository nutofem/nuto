#pragma once

#include "Shape.h"

namespace NuTo
{

class Prism : public Shape
{
public:
    eShape Enum() const override
    {
        return eShape::Prism;
    }

    bool IsWithinShape(Eigen::VectorXd xi, double e = 1e-8) const override
    {
        double x = xi[0];
        double y = xi[1];
        double z = xi[2];
        return (0.-e < x) && (x < 1.+e) && (0.-e < y) && (y < 1.+e) && (-1.-e < z) && (z < 1.+e) && (x + y < 1.+e);
    }
};

} // namespace NuTo
