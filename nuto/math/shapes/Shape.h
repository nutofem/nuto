#pragma once

#include <iosfwd>
#include <Eigen/Core>

namespace NuTo
{

enum class eShape
{
    Line,
    Triangle,
    Quadrilateral,
    Tetrahedron,
    Hexahedron,
    Prism,
    Pyramid
};

class Shape
{
public:
    virtual eShape Enum() const = 0;
    virtual bool IsWithinShape(Eigen::VectorXd xi, double e = 1e-8) const = 0;
    friend std::ostream& operator<<(std::ostream& out, const Shape& shape)
    {
        shape.Info(out);
        return out;
    }
protected:
    virtual void Info(std::ostream& out) const = 0;
};

inline bool operator==(const Shape& lhs, const Shape& rhs)
{
    return lhs.Enum() == rhs.Enum();
}

inline bool operator!=(const Shape& lhs, const Shape& rhs)
{
    return !(lhs == rhs);
}


} // namespace NuTo
