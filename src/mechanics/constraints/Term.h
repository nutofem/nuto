#pragma once
#include "mechanics/nodes/NodeBase.h"

namespace NuTo
{
namespace Constraint
{

class Term
{
public:
    Term(const NodeBase& node, int component, double coefficient)
        : mNode(node)
        , mComponent(component)
        , mCoefficient(coefficient)
    {
    }

    const NodeBase& GetNode() const
    {
        return mNode;
    }

    int GetComponent() const
    {
        return mComponent;
    }

    double GetCoefficient() const
    {
        return mCoefficient;
    }

    void ExchangeNode(const NodeBase& oldNode, const NodeBase& newNode)
    {
        if (&mNode.get() == &oldNode)
            mNode = newNode;
    }

private:
    std::reference_wrapper<const NodeBase> mNode;
    int mComponent;
    double mCoefficient;
};

} /* Constaint */
} /* NuTo */
