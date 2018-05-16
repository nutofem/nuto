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

    bool IsWithinShape(Eigen::VectorXd xi, double e = 1e-8) const override
    {
        double x = xi[0];
        double y = xi[1];
        double z = xi[2];
        return (-(1.-z)-e < x) && (x < (1.-z) + e) && ( -(1.-z)-e < y) && (y< (1.-z) + e);
    }

protected:
    void Info(std::ostream& out) const override
    {
        out << "Pyramid";
    }
};

} // namespace NuTo
