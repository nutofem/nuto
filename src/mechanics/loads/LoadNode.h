#pragma once

#include "mechanics/loads/LoadBase.h"

namespace NuTo
{
class NodeBase;

//! @brief Abstract class for all constraints applied to a single node
class LoadNode : public LoadBase
{

public:
    //! @brief Constructor
    LoadNode(const NodeBase* rNode);

#ifdef ENABLE_SERIALIZATION
    //! @brief serializes (saves) the class
    //! @param ar         archive
    //! @param version    version
    template<class Archive>
    void save(Archive & ar, const unsigned int version) const
    {
    #ifdef DEBUG_SERIALIZATION
        std::cout << "start serialize LoadNode" << std::endl;
    #endif
        ar & boost::serialization::make_nvp("LoadBase", boost::serialization::base_object<LoadBase >(*this));

        std::uintptr_t mNodeAddress = reinterpret_cast<std::uintptr_t>(mNode);
        ar & boost::serialization::make_nvp("mNode", mNodeAddress);
    #ifdef DEBUG_SERIALIZATION
        std::cout << "finish serialize LoadNode" << std::endl;
    #endif
    }

    //! @brief deserializes (loads) the class
    //! @param ar         archive
    //! @param version    version
    template<class Archive>
    void load(Archive & ar, const unsigned int version)
    {
    #ifdef DEBUG_SERIALIZATION
        std::cout << "start serialize LoadNode" << std::endl;
    #endif
        ar & boost::serialization::make_nvp("LoadBase", boost::serialization::base_object<LoadBase >(*this));

        std::uintptr_t mNodeAddress;
        ar & boost::serialization::make_nvp("mNode", mNodeAddress);
        mNode = reinterpret_cast<const NodeBase*>(mNodeAddress);
    #ifdef DEBUG_SERIALIZATION
        std::cout << "finish serialize LoadNode" << std::endl;
    #endif
    }

    BOOST_SERIALIZATION_SPLIT_MEMBER()


    void SetNodePtrAfterSerialization(const std::map<std::uintptr_t, std::uintptr_t>& mNodeMapCast) override
    {
        std::map<std::uintptr_t, std::uintptr_t>::const_iterator it = mNodeMapCast.find(reinterpret_cast<std::uintptr_t>(mNode));
        if (it!=mNodeMapCast.end())
        {
            NodeBase** temp = const_cast<NodeBase**>(&mNode);
            *temp = reinterpret_cast<NodeBase*>(it->second);
        }
        else
            throw Exception("[NuTo::LoadBase::LoadNode] The NodeBase-Pointer could not be updated.");
    }
#endif // ENABLE_SERIALIZATION

protected:
    LoadNode() = default;
    const NodeBase* mNode;
};
} // namespace NuTo
