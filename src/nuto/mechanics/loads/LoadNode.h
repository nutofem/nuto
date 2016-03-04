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
        ar & boost::serialization::make_nvp("LoadNode_LoadBase", boost::serialization::base_object<LoadBase >(*this));
        std::uintptr_t& temp = reinterpret_cast<std::uintptr_t&>(mNode);
        ar & boost::serialization::make_nvp("mNode", temp);
    }

    void SetNodePtrAfterSerialization(const std::map<std::uintptr_t, std::uintptr_t>& mNodeMapCast) override
    {
        std::map<std::uintptr_t, std::uintptr_t>::const_iterator it = mNodeMapCast.find(reinterpret_cast<std::uintptr_t>(mNode));
        if (it!=mNodeMapCast.end())
        {
            NodeBase** temp = const_cast<NodeBase**>(&mNode);
            *temp = reinterpret_cast<NodeBase*>(it->second);
        }
        else
            throw MechanicsException("[NuTo::LoadBase::LoadNode] The NodeBase-Pointer could not be updated.");
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

