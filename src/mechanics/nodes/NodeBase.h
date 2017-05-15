#pragma once

#ifdef ENABLE_SERIALIZATION
#include <boost/serialization/access.hpp>
#include <boost/serialization/export.hpp>
#endif  // ENABLE_SERIALIZATION

#include "mechanics/MechanicsException.h"

#include <eigen3/Eigen/Core>
#include <map>
#include <set>

namespace NuTo
{
namespace Node
{
    enum class eDof : unsigned char;
}// namespace Node

//! @author Thomas Titscher, BAM
//! @date July 2016
//! @brief ... standard abstract class for all nodes
class NodeBase
{
#ifdef ENABLE_SERIALIZATION
    friend class boost::serialization::access;
    //! @brief serializes the class
    //! @param ar         archive
    //! @param version    version
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version);
#endif  // ENABLE_SERIALIZATION

public:
    //! @brief constructor
    NodeBase() {}

    //! @brief destructor
    virtual ~NodeBase() {}

    //! @brief Print information about the node
    friend std::ostream& operator<<(std::ostream& out, const NodeBase& node);

    //! @brief assignment operator
    NodeBase& operator=(NodeBase const& rOther) = default;

    //! @brief sets the global dof number for a specific dof type and a component
    //! @param dof ... specific dof type
    //! @param component ... component index of the dof type (e.g. 0,1,2 for coordinates in 3D)
    //! @param dofNumber ... dof number
    virtual void SetDofNumber(Node::eDof rDof, int rComponent, int dofNumber)
    {
        throw MechanicsException(__PRETTY_FUNCTION__, "Not implemented for this node type.");
    }

    //! @brief returns the number of time derivatives stored at the node
    //! @param rDof ... specific dof type
    //! @return number of derivatives
    virtual int GetNumTimeDerivatives(Node::eDof rDof)const
    {
        throw MechanicsException(__PRETTY_FUNCTION__, "Not implemented for this node type.");
    }

    //! @brief returns if the dof type rDof is an actual degree of freedom (in contrast to COORDINATES)
    virtual bool IsDof(Node::eDof rDof) const
    {
        throw MechanicsException(__PRETTY_FUNCTION__, "Not implemented for this node type.");
    }

    //*************************************************
    //************       ACCESS         ***************
    //*************************************************

    //! @brief returns the total number of dofs
    virtual int GetNumDofs() const
    {
        throw MechanicsException(__PRETTY_FUNCTION__, "Not implemented for this node type.");
    }

    //! @brief returns the number of dofs for a specific dof type
    //! @param rDof ... specific dof type
    virtual int GetNum(Node::eDof rDof) const
    {
        throw MechanicsException(__PRETTY_FUNCTION__, "Not implemented for this node type.");
    }

    //! @brief returns the global dof number for a specific dof type
    //! @param rDof ... specific dof type
    //! @param rComponent ... component index of the dof type (e.g. 0,1,2 for coordinates in 3D)
    virtual int GetDof(Node::eDof rDof, int rComponent) const
    {
        throw MechanicsException(__PRETTY_FUNCTION__, "Not implemented for this node type.");
    }

    //! @brief returns the global dof number for a specific dof type with only one dof (scalar dof)
    //! @param rDof ... specific dof type
    int GetDof(Node::eDof rDof) const;

    //! @brief returns the dof values for a specific dof
    //! @param rDof ... specific dof type
    //! @param rTimeDerivative ... time derivative
    virtual const Eigen::VectorXd& Get(Node::eDof rDof, int rTimeDerivative) const
    {
        throw MechanicsException(__PRETTY_FUNCTION__, "Not implemented for this node type.");
    }

    //! @brief returns the dof values for a specific dof
    //! @param rDof ... specific dof type
    const Eigen::VectorXd& Get(Node::eDof rDof) const
    {
        return Get(rDof, 0);
    }

    //! @brief sets the dof values for a specific dof
    //! @param rDof ... specific dof type
    //! @param rTimeDerivative ... time derivative
    //! @param rValue ... dof value
    virtual void Set(Node::eDof rDof, int rTimeDerivative, const Eigen::VectorXd& rValue)
    {
        throw MechanicsException(__PRETTY_FUNCTION__, "Not implemented for this node type.");
    }

    //! @brief sets the dof values for a specific dof
    //! @param rDof ... specific dof type
    //! @param rValue ... dof value
    void Set(Node::eDof rDof, const Eigen::VectorXd& rValue)
    {
        Set(rDof, 0, rValue);
    }

    //! @brief sets the dof values for a specific scalar dof
    //! @param rDof ... specific dof type
    //! @param rTimeDerivative ... time derivative
    //! @param rValue ... dof scalar value
    void Set(Node::eDof rDof, int rTimeDerivative, double rValue)
    {
        Eigen::VectorXd value(1);
        value[0] = rValue;
        Set(rDof, 0, value);
    }

    //! @brief sets the dof values for a specific scalar dof
    //! @param rDof ... specific dof type
    //! @param rValue ... dof scalar value
    void Set(Node::eDof rDof, double rValue)
    {
        Set(rDof, 0, rValue);
    }

    //! @brief returns a set containing all dof types
    virtual std::set<Node::eDof> GetDofTypes() const
    {
        throw MechanicsException(__PRETTY_FUNCTION__, "Not implemented for this node type.");
    }

    //! @brief clones (copies) the node with all its data, it's supposed to be a new node, so be careful with ptr
    virtual NodeBase* Clone()const=0;

protected:
    //! @brief Outstream function for "virtual friend idiom"
    virtual void Info(std::ostream& out) const {};
};

std::ostream& operator<<(std::ostream& out, const NodeBase& node);

}//namespace NuTo

#ifdef ENABLE_SERIALIZATION
BOOST_CLASS_EXPORT_KEY(NuTo::NodeBase)
#endif // ENABLE_SERIALIZATION
