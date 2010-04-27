#ifndef NODE_GRIDCOORDINATESBASE_H
#define NODE_GRIDCOORDINATESBASE_H

#include "nuto/math/FullMatrix.h"
#include "nuto/mechanics/nodes/NodeBase.h"

namespace NuTo
{
//! @author Andrea Keßler, ISM
//! @date March 2010
//! @brief ... standard class for nodes without coordinates, but a ID (number of the node in the original grid)
class NodeGridCoordinatesBase : public virtual NodeBase
{
#ifdef ENABLE_SERIALIZATION
    friend class boost::serialization::access;
#endif  // ENABLE_SERIALIZATION

public:

    //! @brief constructor
    NodeGridCoordinatesBase() : NodeBase ()
    {
    }

#ifdef ENABLE_SERIALIZATION
    //! @brief serializes the class
    //! @param ar         archive
    //! @param version    version
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version)
    {
    }
#endif // ENABLE_SERIALIZATION

    //! @brief returns the number of coordinates of the node
    //! @return number of coordinates
    virtual int GetNumCoordinates()const=0;

    //! @brief writes the coordinates of a node to the prescribed pointer
    //! @param rCoordinates coordinates
    virtual void GetCoordinates()const=0;

    //! @brief get the ID
     //! @param rNodeNumber number of the node
    virtual int GetNodeId() const=0;


};
}//namespace NuTo
#endif //NODE_GRIDCOORDINATESBASE_H
