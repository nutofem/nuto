#pragma once

#include "Shape.h"

namespace NuTo
{

class Hexahedron : public Shape
{
public:
    eShape Enum() const override
    {
        return eShape::Hexahedron;
    }
};

} // namespace NuTo
