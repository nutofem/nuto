// $Id$

#ifndef NODE_ANGULARVELOCITIES_2D_H_
#define NODE_ANGULARVELOCITIES_2D_H_

#include "nuto/mechanics/nodes/NodeBase.h"

namespace NuTo
{
//! @author Jörg F. Unger
//! @date March 2011
//! @brief ... standard class for nodes having angular velocities degrees of freedom
class NodeAngularVelocities2D : public virtual NuTo::NodeBase
{
#ifdef ENABLE_SERIALIZATION
    friend class boost::serialization::access;
#endif  // ENABLE_SERIALIZATION

public:
    //! @brief default constructor
    NodeAngularVelocities2D();

    //! @brief constructor
    //! @param rVelocities ... prescribed velocity
    NodeAngularVelocities2D(const double rVelocities[1]);

    //! @brief returns the number of angular velocities of the node
    //! @return number of velocities
    int GetNumAngularVelocities()const;

    //! @brief set the angular velocities
    //! @param rVelocities  given velocity
    void SetAngularVelocities2D(const double rVelocities[1]);

    //! @brief writes the velocities of a node to the prescribed pointer
    //! @param rVelocities velocity
    void GetAngularVelocities2D(double rVelocities[1])const;

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
    double mAngularVelocities[1]; //!< node angular velocity
};

}


#endif // NODE_ANGULARVELOCITIES_2D_H_
