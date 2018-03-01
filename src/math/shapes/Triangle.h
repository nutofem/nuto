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
};

} // namespace NuTo
