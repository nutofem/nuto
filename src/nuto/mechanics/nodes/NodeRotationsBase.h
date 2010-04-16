// $Id: $
#ifndef NodeRotationsBase_H
#define NodeRotationsBase_H

#include "nuto/math/FullMatrix.h"
#include "nuto/mechanics/nodes/NodeBase.h"

namespace NuTo
{
//! @author Jörg F. Unger, ISM
//! @date October 2009
//! @brief ... standard class for nodes having rotational degrees of freedom
class NodeRotationsBase : public virtual  NodeBase
{
#ifdef ENABLE_SERIALIZATION
    friend class boost::serialization::access;
#endif // ENABLE_SERIALIZATION
public:

    //! @brief constructor
    NodeRotationsBase() : NodeBase ()
    {
    }

#ifdef ENABLE_SERIALIZATION
    //! @brief serializes the class
    //! @param ar         archive
    //! @param version    version
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version)
    {}
#endif // ENABLE_SERIALIZATION
    //! @brief returns the number of rotations of the node
    //! @return number of rotations
    virtual int GetNumRotations()const=0;

    //! @brief gives the global DOF of a rotation component
    //! @param rComponent component
    //! @return global DOF
    virtual int GetDofRotation(int rComponent)const=0;

    //! @brief set the rotations
    //! @param rRotations  given rotations
    virtual void SetRotations(const double *rRotations)=0;

    //! @brief writes the rotations of a node to the prescribed pointer
    //! @param rRotations rotations
    virtual void GetRotations(double *rRotations)const=0;

protected:
};
}

#endif //NodeRotationsBase_H
