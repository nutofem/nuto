/*
 * BoundaryElementBase.h
 *
 *  Created on: 5 Jun 2015
 *      Author: ttitsche
 */

#ifndef BOUNDARYELEMENTBASE_H_
#define BOUNDARYELEMENTBASE_H_

#include "nuto/mechanics/elements/ElementBase.h"


namespace NuTo
{

namespace BoundaryType
{
enum eType
{
    NOT_SET,
    NEUMANN_HOMOGENEOUS,            // grad nonlocal eq strain * n = 0
    ROBIN_INHOMOGENEOUS,            // l * grad nonlocal eq strain * n + nonlocal eq strain = local eq strain
    MACAULAY                        // l * grad nonlocal eq strain * n + (nonlocal eq strain - local eq strain)_- = 0

};
}

class BoundaryElementBase: public ElementBase
{
public:
    BoundaryElementBase(const ElementBase* rBaseElement, int rSurfaceId);


    virtual ~BoundaryElementBase()
    {

    }

    //! @brief returns the number of nodes in this element
    //! @return number of nodes
    int GetNumNodes() const;

    //! @brief returns a pointer to the i-th node of the element
    //! @param local node number
    //! @return pointer to the node
    NodeBase* GetNode(int rLocalNodeNumber) override;

    //! @brief returns a pointer to the i-th node of the element
    //! @param local node number
    //! @return pointer to the node
    const NodeBase* GetNode(int rLocalNodeNumber) const override;

    //! @brief returns the number of nodes in this element of a specific dof
    //! @brief rDofType dof type
    //! @return number of nodes
    virtual int GetNumNodes(Node::eAttributes rDofType) const override;

    //! @brief returns a pointer to the i-th node of the element
    //! @param local node number
    //! @brief rDofType dof type
    //! @return pointer to the node
    NodeBase* GetNode(int rLocalNodeNumber, Node::eAttributes rDofType) override;

    //! @brief returns a pointer to the i-th node of the element
    //! @param local node number
    //! @brief rDofType dof type
    //! @return pointer to the node
    const NodeBase* GetNode(int rLocalNodeNumber, Node::eAttributes rDofType) const override;

    //! @brief sets the rLocalNodeNumber-th node of the element
    //! @param local node number
    //! @param pointer to the node
    void SetNode(int rLocalNodeNumber, NodeBase* rNode) override
    {
        throw MechanicsException("[NuTo::BoundaryElement1D::SetNode] Probably not needed.");
    }

    //! @brief resizes the node vector
    //! @param rNewNumNodes new number of nodes
    void ResizeNodes(int rNewNumNodes)
    {
        throw MechanicsException("[NuTo::BoundaryElement1D::ResizeNodes] Probably not needed.");
    }

    //! brief exchanges the node ptr in the full data set (elements, groups, loads, constraints etc.)
    //! this routine is used, if e.g. the data type of a node has changed, but the restraints, elements etc. are still identical
    void ExchangeNodePtr(NodeBase* rOldPtr, NodeBase* rNewPtr) override
    {
        throw MechanicsException("[NuTo::BoundaryElement1D::ExchangeNodePtr] Probably not needed.");
    }

    //! @brief sets the section of an element
    //! @param rSection pointer to section
    void SetSection(const SectionBase* rSection) override;

    //! @brief returns a pointer to the section of an element
    //! @return pointer to section
    const SectionBase* GetSection() const override;

    //! @brief calculates the volume of an integration point (weight * detJac)
    //! @return rVolume  vector for storage of the ip volumes (area in 2D, length in 1D)
    const Eigen::VectorXd GetIntegrationPointVolume() const override;

    //! @brief returns the coordinates of an integration point
    //! @param rIpNum integration point
    //! @return rCoordinates coordinates to be returned
    const Eigen::Vector3d GetGlobalIntegrationPointCoordinates(int rIpNum) const override;

    const Eigen::MatrixXd ExtractNodeValues(int rTimeDerivative, Node::eAttributes rAttribute) const override;

protected:
    //! @brief ... just for serialization
    BoundaryElementBase()
    {
    }


    //! @brief ... reorder nodes such that the sign of the length/area/volume of the element changes
    void ReorderNodes() override
    {
        throw MechanicsException("[NuTo::BoundaryElement1D::ReorderNodes] Probably not needed.");
    }

    //! @brief ... check if the element is properly defined (check node dofs, nodes are reordered if the element length/area/volum is negative)
    void CheckElement() override
    {
        throw MechanicsException("[NuTo::BoundaryElement1D::CheckElement] Probably not needed.");
    }

    //The real boundary element that is attached to the virtual boundary element
    const ElementBase* mBaseElement;

    // edge number 0.. left, 1.. right
    int mSurfaceId;

    BoundaryType::eType mBoundaryConditionType;
};

} /* namespace NuTo */

#endif /* BOUNDARYELEMENTBASE_H_ */
