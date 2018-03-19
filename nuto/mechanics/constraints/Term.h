#pragma once
#include "nuto/base/Exception.h"
#include "nuto/mechanics/nodes/NodeSimple.h"

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
    //! @param component component in the dof vector of the node, has to be smaller than node.GetNumValues()
    //! @param coefficient coefficient of the equation term
    Term(const NodeSimple& node, int component, double coefficient)
        : mNode(node)
        , mComponent(component)
        , mCoefficient(coefficient)
    {
        if (component >= node.GetNumValues())
            throw Exception(__PRETTY_FUNCTION__, "Term construction failed. Node has " +
                                                         std::to_string(node.GetNumValues()) +
                                                         " components and you tried to constrain component " +
                                                         std::to_string(component) + ".");
    }

    //! @brief getter for mNode
    const NodeSimple& GetNode() const
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

    int GetConstrainedDofNumber() const
    {
        return mNode.get().GetDofNumber(mComponent);
    }

private:
    //! @brief node reference
    //! @remark `std::reference_wrapper` is used instead of a reference to
    //! make `Term` default CopyAssignable
    std::reference_wrapper<const NodeSimple> mNode;

    //! @brief component component in the dof vector of the node
    int mComponent;

    //! @brief coefficient coefficient of the equation term
    double mCoefficient;
};

} /* Constaint */
} /* NuTo */
