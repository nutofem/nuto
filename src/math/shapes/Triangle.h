#pragma once

#include "Shape.h"
#include "ShapeVisitor.h"

namespace NuTo
{

class Triangle : public Shape
{
public:
    void accept(ShapeVisitor& visitor) const override
    {
        visitor.visit(*this);
    }
};

} // namespace NuTo
