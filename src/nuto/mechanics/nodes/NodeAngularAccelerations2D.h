// $Id$

#ifndef NODE_ANGULAR_ACCELERATIONS_2D_H_
#define NODE_ANGULAR_ACCELERATIONS_2D_H_

#include "nuto/mechanics/nodes/NodeBase.h"

namespace NuTo
{
//! @author Jörg F. Unger
//! @date March 2012
//! @brief ... standard class for nodes having angular accelerations degrees of freedom
class NodeAngularAccelerations2D : public virtual NuTo::NodeBase
{
#ifdef ENABLE_SERIALIZATION
    friend class boost::serialization::access;
#endif  // ENABLE_SERIALIZATION

public:
    //! @brief default constructor
    NodeAngularAccelerations2D();

    //! @brief constructor
    //! @param rAccelerations ... given accelerations
    NodeAngularAccelerations2D(const double rAccelerations[2]);

    //! @brief returns the number of accelerations of the node
    //! @return number of accelerations
    int GetNumAngularAccelerations()const;

    //! @brief set the accelerations
    //! @param rAccelerations  given accelerations
    virtual void SetAngularAccelerations2D(const double rAngularAccelerations[1]);

    //! @brief writes the accelerations of a node to the prescribed pointer
    //! @param rAccelerations accelerations
    void GetAngularAccelerations2D(double rAngularAccelerations[1])const;

    //! @brief returns the type of the node
    //! @return type
    virtual std::string GetNodeTypeStr()const;

    //! @brief returns the type of node as an enum (all the data stored at the node)
    //! @return enum
    virtual Node::eNodeType GetNodeType()const;

#ifdef ENABLE_SERIALIZATION
    //! @brief serializes the class
    //! @param ar         archive
    //! @param version    version
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version);
#endif // ENABLE_SERIALIZATION
protected:
    double mAngularAccelerations[1]; //!< node angular accelerations
};
}
#endif // NODE_ACCELERATIONS_2D_H_
