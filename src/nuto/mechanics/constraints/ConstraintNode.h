// $Id$

#ifndef CONSTRAINTNODE_H
#define CONSTRAINTNODE_H

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

#ifdef ENABLE_SERIALIZATION
    //! @brief serializes the class
    //! @param ar         archive
    //! @param version    version
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version);
#endif // ENABLE_SERIALIZATION

protected:
    //! @brief just for serialization
    ConstraintNode(){};

    const NodeBase* mNode;
};
}//namespace NuTo

#ifdef ENABLE_SERIALIZATION
BOOST_CLASS_EXPORT_KEY(NuTo::ConstraintNode)
#endif // ENABLE_SERIALIZATION

#endif //CONSTRAINTNODE_H

