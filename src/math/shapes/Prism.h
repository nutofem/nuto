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

    bool IsWithinShape(Eigen::VectorXd xi) const override
    {
        double x = xi[0];
        double y = xi[1];
        double z = xi[2];
        return (-1.<x)&&(x<1.)&&(-1.<y)&&(y<1.)&&(0.<z)&&(z<1.)&&(y+z<1.);
    }
};

} // namespace NuTo
