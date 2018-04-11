#pragma once

#include "Shape.h"

namespace NuTo
{

class Spline : public Shape
{
public:
    eShape Enum() const override
    {
        return eShape::Spline;
    }

    bool IsWithinShape(Eigen::VectorXd xi, double e = 1e-8) const override
    {
        (void)e;
        for(int i = 0; i < xi.rows(); i++)
            if(!(0. <= xi(i) && xi(i) <= 1.)) return false;

        return false;
    }
};
} // namespace NuTo
