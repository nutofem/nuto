// $Id$

#pragma once

#ifdef ENABLE_SERIALIZATION
#include <boost/serialization/access.hpp>
#include <boost/serialization/export.hpp>
#include <boost/serialization/split_member.hpp>
#endif  // ENABLE_SERIALIZATION


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
    //! @brief serializes (save) the class
    //! @param ar         archive
    //! @param version    version
    template<class Archive>
    void save(Archive & ar, const unsigned int version) const
    {
    #ifdef DEBUG_SERIALIZATION
        std::cout << "start serialize ConstraintNode" << std::endl;
    #endif
    // The commented code gives warnings in release compile, so split of serialize in save and load
    // needs to be done to avoid those
    /* std::uintptr_t& mNodeAddress = reinterpret_cast<std::uintptr_t&>(mNode);
       ar & boost::serialization::make_nvp("mNode", mNodeAddress);*/

        std::uintptr_t mNodeAddress = reinterpret_cast<std::uintptr_t>(mNode);
        ar & boost::serialization::make_nvp("mNode", mNodeAddress);
    #ifdef DEBUG_SERIALIZATION
        std::cout << "finish serialize ConstraintNode" << std::endl;
    #endif
    }

    //! @brief deserializes (loads) the class
    //! @param ar         archive
    //! @param version    version
    template<class Archive>
    void load(Archive & ar, const unsigned int version)
    {
    #ifdef DEBUG_SERIALIZATION
        std::cout << "start serialize ConstraintNode" << std::endl;
    #endif
        std::uintptr_t mNodeAddress;
        ar & boost::serialization::make_nvp("mNode", mNodeAddress);
        mNode = reinterpret_cast<const NodeBase*>(mNodeAddress);
    #ifdef DEBUG_SERIALIZATION
        std::cout << "finish serialize ConstraintNode" << std::endl;
    #endif
    }

    BOOST_SERIALIZATION_SPLIT_MEMBER()

    //! @brief NodeBase-Pointer are not serialized to avoid cyclic dependencies, but are serialized as Pointer-Address (uintptr_t)
    //! Deserialization of the NodeBase-Pointer is done by searching and casting back the address in the map
    //! @param mNodeMapCast std::map containing the old and new addresses
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

