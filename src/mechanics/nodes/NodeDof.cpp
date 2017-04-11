#include "mechanics/nodes/NodeDof.h"
#include "mechanics/nodes/NodeEnum.h"

using namespace NuTo;

#ifdef ENABLE_SERIALIZATION
#include "math/CustomBoostSerializationExtensions.h"
#include <boost/serialization/unordered_map.hpp>
#include <boost/serialization/vector.hpp>
#endif // ENABLE_SERIALIZATION

NodeDof::NodeDof(std::map<Node::eDof, NodeDofInfo> rDofInfos)
    : NodeBase::NodeBase()
{
    for (auto& it : rDofInfos)
    {
        Node::eDof dofType = it.first;
        const NodeDofInfo& info = it.second;

        // allocate global dof numbers
        if (info.mIsDof)
        {
            mDofNumbers[dofType] = Eigen::VectorXi::Zero(info.mDimension);
        }

        // allocate dof values
        int vectorSize = info.mNumTimeDerivatives + 1; // +1 since 0th time derivative --> size 1, ...
        mDofValues[dofType].resize(vectorSize);
        for (int i = 0; i < vectorSize; ++i)
            mDofValues[dofType][i] = Eigen::VectorXd::Zero(info.mDimension);
    }
}


void NodeDof::SetDofNumber(Node::eDof dof, int component, int dofNumber)
{
    mDofNumbers[dof][component] = dofNumber;
}

int NodeDof::GetNumTimeDerivatives(Node::eDof rDof) const
{
    const auto& it = mDofValues.find(rDof);
    if (it == mDofValues.end())
        return 0;

    return it->second.size() - 1; // -1: dt 0 --> size 1; dt 1 --> size 2; ...
}

bool NodeDof::IsDof(Node::eDof rDof) const
{
    return mDofNumbers.find(rDof) != mDofNumbers.end();
}


int NodeDof::GetNumDofs() const
{
    int numDofs = 0;
    for (auto& it : mDofNumbers)
        numDofs += it.second.rows();
    return numDofs;
}


int NodeDof::GetNum(Node::eDof rDof) const
{
    const auto& it = mDofValues.find(rDof);
    if (it == mDofValues.end())
        return 0;

    return it->second[0].rows();
}

int NodeDof::GetDof(Node::eDof rDof, int rComponent) const
{
    const auto& it = mDofNumbers.find(rDof);
    if (it == mDofNumbers.end())
        throw MechanicsException(__PRETTY_FUNCTION__, "Cannot access dof type " + Node::DofToString(rDof));

    if (rComponent >= it->second.size())
        throw MechanicsException(__PRETTY_FUNCTION__, "Cannot access component " + std::to_string(rComponent) +
                                                              ". Component " + std::to_string(it->second.size()) +
                                                              " was requested.");

    return it->second[rComponent];
}

const Eigen::VectorXd& NodeDof::Get(Node::eDof rDof, int rTimeDerivative) const
{
    const auto& it = mDofValues.find(rDof);
    if (it == mDofValues.end())
        throw MechanicsException(__PRETTY_FUNCTION__, "Cannot access dof type " + Node::DofToString(rDof));

    if (rTimeDerivative >= (int)it->second.size())
        throw MechanicsException(__PRETTY_FUNCTION__,
                                 "Cannot access time derivative " + std::to_string(rTimeDerivative) +
                                         ". This node only has " + std::to_string(it->second.size()) + ".");

    return it->second[rTimeDerivative];
}

void NodeDof::Set(Node::eDof rDof, int rTimeDerivative, const Eigen::VectorXd& rValue)
{
    auto it = mDofValues.find(rDof);
    if (it == mDofValues.end())
        throw MechanicsException(__PRETTY_FUNCTION__, "Cannot access dof type " + Node::DofToString(rDof));

    if (rTimeDerivative >= (int)it->second.size())
        throw MechanicsException(__PRETTY_FUNCTION__,
                                 "Cannot access time derivative " + std::to_string(rTimeDerivative) +
                                         ". This node only has " + std::to_string(it->second.size()) + ".");

    it->second[rTimeDerivative] = rValue;
}

std::set<Node::eDof> NodeDof::GetDofTypes() const
{
    std::set<Node::eDof> dofTypes;
    for (auto it : mDofValues)
        dofTypes.insert(it.first);
    return dofTypes;
}

void NodeDof::Info(std::ostream& out) const
{
    for (auto& it : mDofValues)
    {
        out << Node::DofToString(it.first) << ": " << it.second[0].rows() << " dt:" << it.second.size() << "\n";
    }
}

NodeBase* NodeDof::Clone() const
{
    return new NodeDof(*this);
}

#ifdef ENABLE_SERIALIZATION
//! @brief serializes the class
//! @param ar         archive
//! @param version    version
template void NodeDof::serialize(boost::archive::binary_oarchive& ar, const unsigned int version);
template void NodeDof::serialize(boost::archive::binary_iarchive& ar, const unsigned int version);
template void NodeDof::serialize(boost::archive::xml_oarchive& ar, const unsigned int version);
template void NodeDof::serialize(boost::archive::xml_iarchive& ar, const unsigned int version);
template void NodeDof::serialize(boost::archive::text_oarchive& ar, const unsigned int version);
template void NodeDof::serialize(boost::archive::text_iarchive& ar, const unsigned int version);
template <class Archive>
void NodeDof::serialize(Archive& ar, const unsigned int version)
{
#ifdef DEBUG_SERIALIZATION
    std::cout << "start serialize NodeDof"
              << "\n";
#endif
    ar& BOOST_SERIALIZATION_BASE_OBJECT_NVP(NodeBase) & BOOST_SERIALIZATION_NVP(mDofValues) &
            BOOST_SERIALIZATION_NVP(mDofNumbers);
#ifdef DEBUG_SERIALIZATION
    std::cout << "finish serialize NodeDof \n";
#endif
}
BOOST_CLASS_EXPORT_IMPLEMENT(NodeDof)
#endif // ENABLE_SERIALIZATION
