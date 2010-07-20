// $Id$

#ifndef NODE_VELOCITIES_3D_H_
#define NODE_VELOCITIES_3D_H_

#include "nuto/mechanics/nodes/NodeBase.h"

namespace NuTo
{
//! @author Stefan Eckardt, IFF
//! @date July 2010
//! @brief ... standard class for nodes having velocities degrees of freedom
class NodeVelocities3D : public virtual NuTo::NodeBase
{
#ifdef ENABLE_SERIALIZATION
    friend class boost::serialization::access;
#endif  // ENABLE_SERIALIZATION

public:
    //! @brief default constructor
    NodeVelocities3D();

    //! @brief constructor
    //! @param rVelocities ... prescribed velocity
    NodeVelocities3D(const double rVelocities[3]);

    //! @brief returns the number of velocities of the node
    //! @return number of velocities
    int GetNumVelocities()const;

    //! @brief set the velocities
    //! @param rVelocities  given velocity
    void SetVelocities3D(const double rVelocities[3]);

    //! @brief writes the velocities of a node to the prescribed pointer
    //! @param rVelocities velocity
    void GetVelocities3D(double rVelocities[3])const;

    //! @brief returns the type of the node
    //! @return type
    std::string GetNodeTypeStr()const;

#ifdef ENABLE_SERIALIZATION
    //! @brief serializes the class
    //! @param ar         archive
    //! @param version    version
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version);
#endif // ENABLE_SERIALIZATION
protected:
    double mVelocities[3]; //!< node velocity
};
}
#endif // NODE_VELOCITIES_3D_H_
