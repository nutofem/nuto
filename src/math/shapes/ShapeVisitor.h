#pragma once

namespace NuTo
{

class Triangle;
class Hexahedron;

class ShapeVisitor
{
public:
    virtual void visit(const Triangle&) = 0;
    virtual void visit(const Hexahedron&) = 0;
};

} // namespace NuTo
