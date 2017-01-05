#pragma once

#ifdef ENABLE_SERIALIZATION
#include <boost/serialization/access.hpp>
#include <boost/serialization/export.hpp>
#endif  // ENABLE_SERIALIZATION

#ifdef ENABLE_VISUALIZE
#include <list>
#include <memory>
#endif // ENABLE_VISUALIZE



#include "mechanics/MechanicsException.h"

#include <eigen3/Eigen/Core>
#include <map>
#include <set>
#include <vector>


namespace NuTo
{
#ifdef ENABLE_VISUALIZE
class VisualizeComponent;
class VisualizeUnstructuredGrid;
#endif // ENABLE_VISUALIZE

template <class T, int rows> class FullVector;

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

    //! @brief assignment operator
    NodeBase& operator=(NodeBase const& rOther) = default;

    //! @brief sets the global dofs numbers for each dof type
    //! @param rDofNumbers ... map containing the dof type and the current number
    virtual void SetGlobalDofsNumbers(std::map<Node::eDof, int>& rDofNumbers)
    {
        throw MechanicsException(__PRETTY_FUNCTION__, "Not implemented for node type " + GetNodeTypeStr() + ".");
    }

    //! @brief write dof values to the node (based on global dof number) for a specific dof type
    //! @param rTimeDerivative ... time derivative (e.g. 0 disp, 1 vel, 2 acc)
    //! @param rDofType ... specific dof type
    //! @param rActiveDofValues ... active dof values
    //! @param rDependentDofValues ... dependent dof values
    virtual void SetGlobalDofValues(int rTimeDerivative, Node::eDof rDofType, const Eigen::VectorXd& rActiveDofValues, const Eigen::VectorXd& rDependentDofValues)
    {
        throw MechanicsException(__PRETTY_FUNCTION__, "Not implemented for node type " + GetNodeTypeStr() + ".");
    }

    //! @brief extract dof values from the node (based on global dof number) for a specific dof type
    //! @param rTimeDerivative ... time derivative (e.g. 0 disp, 1 vel, 2 acc)
    //! @param rDofType ... specific dof type
    //! @param rActiveDofValues ... active dof values
    //! @param rDependentDofValues ... dependent dof values
    virtual void GetGlobalDofValues(int rTimeDerivative, Node::eDof rDofType, Eigen::VectorXd& rActiveDofValues, Eigen::VectorXd& rDependentDofValues) const
    {
        throw MechanicsException(__PRETTY_FUNCTION__, "Not implemented for node type " + GetNodeTypeStr() + ".");
    }

    //! @brief renumber the global dofs according to predefined ordering
    //! @param rMappingInitialToNewOrdering ... mapping from initial ordering to the new ordering
    virtual void RenumberGlobalDofs(Node::eDof rDof, std::vector<int>& rMappingInitialToNewOrdering)
    {
        throw MechanicsException(__PRETTY_FUNCTION__, "Not implemented for node type " + GetNodeTypeStr() + ".");
    }


    //! @brief returns the number of time derivatives stored at the node
    //! @param rDof ... specific dof type
    //! @return number of derivatives
    virtual int GetNumTimeDerivatives(Node::eDof rDof)const
    {
        throw MechanicsException(__PRETTY_FUNCTION__, "Not implemented for node type " + GetNodeTypeStr() + ".");
    }

    //! @brief returns if the dof type rDof is an actual degree of freedom (in contrast to COORDINATES)
    virtual bool IsDof(Node::eDof rDof) const
    {
        throw MechanicsException(__PRETTY_FUNCTION__, "Not implemented for node type " + GetNodeTypeStr() + ".");
    }

    //*************************************************
    //************       ACCESS         ***************
    //*************************************************

    //! @brief returns the total number of dofs
    virtual int GetNumDofs() const
    {
        throw MechanicsException(__PRETTY_FUNCTION__, "Not implemented for node type " + GetNodeTypeStr() + ".");
    }

    //! @brief returns the number of dofs for a specific dof type
    //! @param rDof ... specific dof type
    virtual int GetNum(Node::eDof rDof) const
    {
        throw MechanicsException(__PRETTY_FUNCTION__, "Not implemented for node type " + GetNodeTypeStr() + ".");
    }

    //! @brief returns the global dof number for a specific dof type
    //! @param rDof ... specific dof type
    //! @param rComponent ... component index of the dof type (e.g. 0,1,2 for coordinates in 3D)
    virtual int GetDof(Node::eDof rDof, int rComponent) const
    {
        throw MechanicsException(__PRETTY_FUNCTION__, "Not implemented for node type " + GetNodeTypeStr() + ".");
    }

    //! @brief returns the global dof number for a specific dof type with only one dof (scalar dof)
    //! @param rDof ... specific dof type
    int GetDof(Node::eDof rDof) const;

    //! @brief returns the dof values for a specific dof
    //! @param rDof ... specific dof type
    //! @param rTimeDerivative ... time derivative
    virtual const Eigen::VectorXd& Get(Node::eDof rDof, int rTimeDerivative) const
    {
        throw MechanicsException(__PRETTY_FUNCTION__, "Not implemented for node type " + GetNodeTypeStr() + ".");
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
        throw MechanicsException(__PRETTY_FUNCTION__, "Not implemented for node type " + GetNodeTypeStr() + ".");
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
        throw MechanicsException(__PRETTY_FUNCTION__, "Not implemented for node type " + GetNodeTypeStr() + ".");
    }


    //! @brief returns the type of node as a string (all the data stored at the node)
    //! @return string
    virtual std::string GetNodeTypeStr()const=0;



    //! @brief clones (copies) the node with all its data, it's supposed to be a new node, so be careful with ptr
    virtual NodeBase* Clone()const=0;

#ifdef ENABLE_VISUALIZE
    virtual void Visualize(VisualizeUnstructuredGrid& rVisualize, const std::list<std::shared_ptr<NuTo::VisualizeComponent>>& rVisualizationList) const;


private:
    //! @brief extracts the desired vector data to a 3d vector
    Eigen::Vector3d GetPointVectorData(Node::eDof rDofType, int rTimeDerivative) const;
#endif // ENABLE_VISUALIZE


protected:
    //the base class of the nodes must not contain any data



};

class less_XCoordinate2D : public std::binary_function<NodeBase*, NodeBase* , bool>
{
public:

    //! @brief sorts the nodes in increasing x-direction
    less_XCoordinate2D()
    {
    }

    bool operator()(NodeBase* nodePtr1, NodeBase* nodePtr2);
};

class greater_XCoordinate2D : public std::binary_function<NodeBase*, NodeBase* , bool>
{
public:
    greater_XCoordinate2D()
    {
    }

    bool operator()(NodeBase* nodePtr1, NodeBase* nodePtr2);
};

class less_YCoordinate2D : public std::binary_function<NodeBase*, NodeBase* , bool>
{
public:
    less_YCoordinate2D()
    {
    }

    bool operator()(NodeBase* nodePtr1, NodeBase* nodePtr2);
};

class greater_YCoordinate2D : public std::binary_function<NodeBase*, NodeBase* , bool>
{
public:
    greater_YCoordinate2D()
    {
    }

    bool operator()(NodeBase* nodePtr1, NodeBase* nodePtr2);
};

}//namespace NuTo

#ifdef ENABLE_SERIALIZATION
BOOST_CLASS_EXPORT_KEY(NuTo::NodeBase)
#endif // ENABLE_SERIALIZATION
