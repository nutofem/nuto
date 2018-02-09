#pragma once

namespace NuTo
{

class ShapeVisitor;

class Shape
{
public:
    virtual void accept(ShapeVisitor& visitor) const = 0;
};

} // namespace NuTo
