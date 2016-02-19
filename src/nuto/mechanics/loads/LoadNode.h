// $Id$
#ifndef LOADNODE_H
#define LOADNODE_H

#ifdef ENABLE_SERIALIZATION
#include <boost/serialization/access.hpp>
#include <boost/serialization/export.hpp>
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/xml_oarchive.hpp>
#include <boost/archive/xml_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#endif  // ENABLE_SERIALIZATION

#include "nuto/mechanics/MechanicsException.h"
#include "nuto/mechanics/loads/LoadBase.h"

namespace NuTo
{
class NodeBase;
//! @author Jörg F. Unger, ISM
//! @date October 2009
//! @brief ... abstract class for all constraints applied to a single node
class LoadNode : public LoadBase
{
#ifdef ENABLE_SERIALIZATION
    friend class boost::serialization::access;
#endif  // ENABLE_SERIALIZATION

public:
    //! @brief constructor
    LoadNode(){}

    //! @brief constructor
    LoadNode(int rLoadCase, const NodeBase* rNode);

#ifdef ENABLE_SERIALIZATION
    //! @brief serializes the class
    //! @param ar         archive
    //! @param version    version
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version)
    {
//        ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(LoadBase);
        ar & boost::serialization::make_nvp("LoadNode_LoadBase", boost::serialization::base_object<LoadBase >(*this));
//        ar & BOOST_SERIALIZATION_NVP(const_cast<NodeBase*&>(mNode));
        ar & boost::serialization::make_nvp ("LoadNode_mNode", const_cast<NodeBase*&>(mNode));
    }
#endif // ENABLE_SERIALIZATION

protected:
    const NodeBase* mNode;
};
}//namespace NuTo

#ifdef ENABLE_SERIALIZATION
BOOST_CLASS_EXPORT_KEY(NuTo::LoadNode)
#endif

#endif //LOADNODE_H

