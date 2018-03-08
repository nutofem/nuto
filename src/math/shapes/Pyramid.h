#pragma once

#include "Shape.h"

namespace NuTo
{

class Pyramid : public Shape
{
public:
    eShape Enum() const override
    {
        return eShape::Pyramid;
    }

    bool IsWithinShape(Eigen::VectorXd xi) const override
    {
        double x = xi[0];
        double y = xi[1];
        double z = xi[2];
        return (-(1.-z)<x)&&(x<(1.-z))&&(-(1.-z)<y)&&(y<(1.-z));
    }
};

} // namespace NuTo
