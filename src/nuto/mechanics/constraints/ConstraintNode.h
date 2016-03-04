// $Id$

#ifndef CONSTRAINTNODE_H
#define CONSTRAINTNODE_H

#ifdef ENABLE_SERIALIZATION
#include <boost/serialization/access.hpp>
#include <boost/serialization/export.hpp>
#endif  // ENABLE_SERIALIZATION

#include "nuto/mechanics/constraints/ConstraintEnum.h"

namespace NuTo
{
class NodeBase;
//! @author Jörg F. Unger, ISM
//! @date October 2009
//! @brief ... abstract class for all constraints applied to a single node
class ConstraintNode
{
#ifdef ENABLE_SERIALIZATION
    friend class boost::serialization::access;
#endif  // ENABLE_SERIALIZATION

public:
    //! @brief constructor
    ConstraintNode(const NodeBase* rNode);

    //! @brief destructor
    virtual ~ConstraintNode();

#ifdef ENABLE_SERIALIZATION
    //! @brief serializes the class
    //! @param ar         archive
    //! @param version    version
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version);

    //! @brief NodeBase-Pointer are not serialized to avoid cyclic dependencies, but are serialized as Pointer-Adress (uintptr_t)
    //! Deserialization of the NodeBase-Pointer is done by searching and casting back the adress in the map
    //! @param mNodeMapCast std::map containing the old and new adresses
    virtual void SetNodePtrAfterSerialization(const std::map<uintptr_t, uintptr_t>& mNodeMapCast);
#endif // ENABLE_SERIALIZATION

protected:
    //! @brief just for serialization
    ConstraintNode(){}

    const NodeBase* mNode;
};
}//namespace NuTo

#ifdef ENABLE_SERIALIZATION
BOOST_CLASS_EXPORT_KEY(NuTo::ConstraintNode)
#endif // ENABLE_SERIALIZATION

#endif //CONSTRAINTNODE_H

