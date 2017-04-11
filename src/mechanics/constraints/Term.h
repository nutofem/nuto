#pragma once
#include "mechanics/nodes/NodeBase.h"

namespace NuTo
{
namespace Constraint
{
//! @brief stores constraint equation terms
class Term
{
public:
    //! @brief ctor that sets all members so that the term consists of
    //! ... + coefficient * node->GetDof(component) + ...
    //! @param node node reference
    //! @param component component in the dof vector of the node
    //! @param coefficient coefficient of the equation term
    Term(const NodeBase& node, int component, double coefficient)
        : mNode(node)
        , mComponent(component)
        , mCoefficient(coefficient)
    {
    }

    //! @brief replaces oldNode with newNode
    //! @param oldNode old node
    //! @param newNode new node
    void ExchangeNode(const NodeBase& oldNode, const NodeBase& newNode)
    {
        if (&mNode.get() == &oldNode)
            mNode = newNode;
    }


    //! @brief getter for mNode
    const NodeBase& GetNode() const
    {
        return mNode;
    }

    //! @brief getter for mComponent
    int GetComponent() const
    {
        return mComponent;
    }

    //! @brief getter for mCoefficient
    double GetCoefficient() const
    {
        return mCoefficient;
    }

private:
    //! @brief node reference, std::reference_wrapper is neat since it provides copy
    //! ctor and assignment
    std::reference_wrapper<const NodeBase> mNode;

    //! @brief component component in the dof vector of the node
    int mComponent;

    //! @brief coefficient coefficient of the equation term
    double mCoefficient;
};

} /* Constaint */
} /* NuTo */
