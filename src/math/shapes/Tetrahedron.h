#pragma once

#include "Shape.h"

namespace NuTo
{

class Tetrahedron : public Shape
{
public:
    eShape Enum() const override
    {
        return eShape::Tetrahedron;
    }
};

} // namespace NuTo
