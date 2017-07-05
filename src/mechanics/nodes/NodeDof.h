#pragma once

#include <vector>
#include "mechanics/nodes/NodeBase.h"
#include "mechanics/nodes/DofHash.h"

namespace NuTo
{
struct NodeDofInfo
{
    int mDimension;
    int mNumTimeDerivatives;
    bool mIsDof;
};


//! @author Thomas Titscher, BAM
//! @date July 2016
//! @brief ... standard class for all nodes
class NodeDof: public NodeBase
{
public:
    //! @brief ctor
    //! @param contains the information for each dof type
    NodeDof(std::map<Node::eDof, NodeDofInfo> rDofInfos);

    //! @brief sets the global dof number for a specific dof type and a component
    //! @param dof ... specific dof type
    //! @param component ... component index of the dof type (e.g. 0,1,2 for coordinates in 3D)
    //! @param dofNumber ... dof number
    void SetDofNumber(Node::eDof dof, int domponent, int dofNumber) override;

    //! @brief returns the number of time derivatives stored at the node
    //! @param rDof ... specific dof type
    //! @return number of derivatives
    int GetNumTimeDerivatives(Node::eDof rDof) const override;

    //! @brief returns if the dof type rDof is an actual degree of freedom (in contrast to COORDINATES)
    bool IsDof(Node::eDof rDof) const override;

    //*************************************************
    //************       ACCESS         ***************
    //*************************************************

    //! @brief returns the total number of dofs
    int GetNumDofs() const override;

    //! @brief returns the number of dofs for a specific dof type
    //! @param rDof ... specific dof type
    int GetNum(Node::eDof rDof) const override;

    //! @brief returns the global dof number for a specific dof type
    //! @param rDof ... specific dof type
    //! @param rComponent ... component index of the dof type (e.g. 0,1,2 for coordinates in 3D)
    int GetDof(Node::eDof rDof, int rComponent) const override;

    //! @brief returns the dof values for a specific dof
    //! @param rDof ... specific dof type
    //! @param rTimeDerivative ... time derivative
    const Eigen::VectorXd& Get(Node::eDof rDof, int rTimeDerivative) const override;

    //! @brief sets the dof values for a specific dof
    //! @param rDof ... specific dof type
    //! @param rTimeDerivative ... time derivative
    //! @param rValue ... dof value
    void Set(Node::eDof rDof, int rTimeDerivative , const Eigen::VectorXd& rValue) override;

    //! @brief returns a set containing all dof types
    std::set<Node::eDof> GetDofTypes() const override;

    //! @brief clones (copies) the node with all its data, it's supposed to be a new node, so be careful with ptr
    NodeBase* Clone() const override;

protected:
    //! @brief Outstream function for "virtual friend idiom"
    virtual void Info(std::ostream& out) const override;

private:

    //! @brief stores the dof values (std::vector for time derivatives, VectorXd for values)
    std::map<Node::eDof, std::vector<Eigen::VectorXd>> mDofValues;

    //! @brief stores the global dof numbers in a dynamic integer vector
    std::map<Node::eDof, Eigen::VectorXi> mDofNumbers;

};

} /* namespace NuTo */

