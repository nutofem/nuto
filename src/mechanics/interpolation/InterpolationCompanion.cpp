#include "InterpolationCompanion.h"

#include "base/Exception.h"

#include "math/shapes/ShapeVisitor.h"
#include "math/shapes/Triangle.h"
#include "math/shapes/Hexahedron.h"

#include "InterpolationTriangleLinear.h"
#include "InterpolationTriangleQuadratic.h"
#include "InterpolationBrickLinear.h"

using namespace NuTo;

class InterpolationFactory : public ShapeVisitor
{
public:
    void visit(const Triangle&) override
    {
        switch (mOrder)
        {
        case 1:
            mOutput = std::make_unique<InterpolationTriangleLinear>();
            break;
        case 2:
            mOutput = std::make_unique<InterpolationTriangleQuadratic>();
            break;
        default:
            throw Exception(__PRETTY_FUNCTION__,
                            "Triangle interpolation of order " + std::to_string(mOrder) + " is not defined.");
        }
    }

    void visit(const Hexahedron&) override
    {
        switch (mOrder)
        {
        case 1:
            mOutput = std::make_unique<InterpolationBrickLinear>();
            break;
        default:
            throw Exception(__PRETTY_FUNCTION__,
                            "Triangle interpolation of order " + std::to_string(mOrder) + " is not defined.");
        }
    }

    std::unique_ptr<InterpolationSimple> Create(const Shape& shape, int order)
    {
        mOrder = order;
        shape.accept(*this);
        return std::move(mOutput);
    }

private:
    std::unique_ptr<InterpolationSimple> mOutput;
    int mOrder;
};

std::unique_ptr<InterpolationSimple> NuTo::CreateInterpolation(const Shape& shape, int order)
{
    InterpolationFactory tmpFactory;
    return tmpFactory.Create(shape, order);
}
