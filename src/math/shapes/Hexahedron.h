#pragma once

#include "Shape.h"
#include "ShapeVisitor.h"

namespace NuTo
{

class Hexahedron : public Shape
{
public:
    void accept(ShapeVisitor& visitor) const override
    {
        visitor.visit(*this);
    }
};

} // namespace NuTo
