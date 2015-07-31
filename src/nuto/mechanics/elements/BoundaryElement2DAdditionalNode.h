#ifndef BOUNDARYELEMENT2DADDITIONALNODE_H
#define BOUNDARYELEMENT2DADDITIONALNODE_H

#include "nuto/mechanics/nodes/NodeBase.h"
#include "nuto/mechanics/elements/BoundaryElement2D.h"
namespace NuTo
{

class BoundaryElement2DAdditionalNode : public BoundaryElement2D
{
public:
    BoundaryElement2DAdditionalNode(const ElementBase *rBaseElement,
                                    int rSurfaceEdge,
                                    NodeBase* rNodeDependent);
/*
    //! @brief calculates output data for the element
    //! @param eOutput ... coefficient matrix 0 1 or 2  (mass, damping and stiffness) and internal force (which includes inertia terms)
    //! @param updateStaticData (with DummyOutput), IPData, globalrow/column dofs etc.
    Error::eError Evaluate(boost::ptr_multimap<NuTo::Element::eOutput, NuTo::ElementOutputBase>& rElementOutput) override;
*/

    //! @brief returns the enum (type of the element)
    //! @return enum
    NuTo::Element::eElementType GetEnumType() const override
    {
        return Element::BOUNDARYELEMENT2DADDITIONALNODE;
    }

    //! @brief Gets the additional node of an boundary element, if it has one
    //! @return Additional boundary node
    virtual NodeBase *GetAdditionalBoundaryNode() const override;

    /*
    //! @brief returns the number of nodes in this element that are influenced by it
    //! @remark overridden by boundary elements
    //! @return number of nodes
    int GetNumInfluenceNodes()const override
    {
        return mBaseElement->GetNumNodes() + 1;
    }

    //! @brief returns a pointer to the i-th node of the element
    //! @remark overridden by boundary elements
    //! @param local node number
    //! @return pointer to the node
    const NodeBase* GetInfluenceNode(int rLocalNodeNumber)const override
    {
        return mBaseElement->GetNode(rLocalNodeNumber);
    }
    */

private:
    NodeBase* mNodeDependent;

};

} //namespace NuTo

#endif // BOUNDARYELEMENT2DADDITIONALNODE_H
