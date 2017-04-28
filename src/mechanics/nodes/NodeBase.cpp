#include "mechanics/nodes/NodeBase.h"
#include "mechanics/nodes/NodeEnum.h"
#include <cassert>

using namespace NuTo;

int NodeBase::GetDof(Node::eDof rDof) const
{
    assert(GetNum(rDof) == 1);
    return GetDof(rDof, 0);
}

namespace NuTo
{
std::ostream& operator<<(std::ostream& out, const NodeBase& node)
{
    node.Info(out);
    return out;
}
} /* NuTo */

#ifdef ENABLE_SERIALIZATION
//! @brief serializes the class
//! @param ar         archive
//! @param version    version
template void NodeBase::serialize(boost::archive::binary_oarchive& ar, const unsigned int version);
template void NodeBase::serialize(boost::archive::xml_oarchive& ar, const unsigned int version);
template void NodeBase::serialize(boost::archive::text_oarchive& ar, const unsigned int version);
template void NodeBase::serialize(boost::archive::binary_iarchive& ar, const unsigned int version);
template void NodeBase::serialize(boost::archive::xml_iarchive& ar, const unsigned int version);
template void NodeBase::serialize(boost::archive::text_iarchive& ar, const unsigned int version);
template <class Archive>
void NodeBase::serialize(Archive& ar, const unsigned int version)
{
#ifdef DEBUG_SERIALIZATION
    std::cout << "start serialize NodeBase"
              << "\n";
#endif
// nothing to do here, no members...
#ifdef DEBUG_SERIALIZATION
    std::cout << "finish serialize NodeBase \n";
#endif
}
BOOST_CLASS_EXPORT_IMPLEMENT(NodeBase)
#endif // ENABLE_SERIALIZATION

