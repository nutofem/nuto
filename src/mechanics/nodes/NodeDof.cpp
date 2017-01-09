

#include "mechanics/nodes/NodeDof.h"
#include "mechanics/nodes/NodeEnum.h"

#ifdef ENABLE_SERIALIZATION
#include "math/CustomBoostSerializationExtensions.h"
#include <boost/serialization/unordered_map.hpp>
#include <boost/serialization/vector.hpp>
#endif // ENABLE_SERIALIZATION

NuTo::NodeDof::NodeDof(std::map<Node::eDof, NodeDofInfo> rDofInfos)
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
        int vectorSize = info.mNumTimeDerivatives +1; // +1 since 0th time derivative --> size 1, ...
        mDofValues[dofType].resize(vectorSize);
        for (int i = 0; i < vectorSize; ++i)
            mDofValues[dofType][i] = Eigen::VectorXd::Zero(info.mDimension);
    }
}



void NuTo::NodeDof::SetGlobalDofsNumbers(std::map<Node::eDof, int>& rDofNumbers)
{
    for (auto& dofTypeNumberPair : mDofNumbers)
    {
        Node::eDof dofType = dofTypeNumberPair.first;
        auto& dofNumbers = dofTypeNumberPair.second;

        for (int i = 0; i < dofNumbers.rows(); ++i)
        {
            dofNumbers[i] = rDofNumbers[dofType]++;
        }
    }
}


void NuTo::NodeDof::SetGlobalDofValues(
        int rTimeDerivative,
        Node::eDof rDofType,
        const Eigen::VectorXd& rActiveDofValues,
        const Eigen::VectorXd& rDependentDofValues)
{
    auto it = mDofValues.find(rDofType);

    if (mDofValues.find(rDofType) == mDofValues.end())
        return; // the node does not have the requested dof type

    assert(GetNumTimeDerivatives(rDofType) >= rTimeDerivative);

    auto& values = it->second[rTimeDerivative];

    for (int i = 0; i < values.rows(); ++i)
    {
        double dofValueToSet = GetDofValueFromVector(mDofNumbers[rDofType][i], rActiveDofValues, rDependentDofValues);
        values[i] = dofValueToSet;
    }
}

void NuTo::NodeDof::GetGlobalDofValues(
        int rTimeDerivative,
        Node::eDof rDofType,
        Eigen::VectorXd& rActiveDofValues,
        Eigen::VectorXd& rDependentDofValues) const
{
    const auto& it = mDofValues.find(rDofType);

    if (it == mDofValues.end())
        return; // the node does not have the requested dof type

    assert(GetNumTimeDerivatives(rDofType) >= rTimeDerivative);

    auto& values = it->second[rTimeDerivative];

    for (int i = 0; i < values.rows(); ++i)
    {
        double dofValueToWrite = values[i];
        WriteNodeValueToVector(mDofNumbers.at(rDofType)[i], dofValueToWrite, rActiveDofValues, rDependentDofValues);
    }
}

void NuTo::NodeDof::RenumberGlobalDofs(Node::eDof rDofType, std::vector<int>& rMappingInitialToNewOrdering)
{
    auto it = mDofNumbers.find(rDofType);
    if(it == mDofNumbers.end())
        return; // This is not an error. E.g. linear disps, quadratic temp. Some nodes only have Temp. Still, disps are called for renumbering.

    auto& dofNumbers = it->second;
    for (int i = 0; i < dofNumbers.rows(); ++i)
    {
        dofNumbers[i] = rMappingInitialToNewOrdering[dofNumbers[i]];
    }
}

int NuTo::NodeDof::GetNumTimeDerivatives(Node::eDof rDof) const
{
    const auto& it = mDofValues.find(rDof);
    if (it == mDofValues.end())
        return 0;

    return it->second.size() - 1; // -1: dt 0 --> size 1; dt 1 --> size 2; ...
}

bool NuTo::NodeDof::IsDof(Node::eDof rDof) const
{
    return mDofNumbers.find(rDof) != mDofNumbers.end();
}


int NuTo::NodeDof::GetNumDofs() const
{
    int numDofs = 0;
    for (auto& it : mDofNumbers)
        numDofs += it.second.rows();
    return numDofs;
}


int NuTo::NodeDof::GetNum(Node::eDof rDof) const
{
    const auto& it = mDofValues.find(rDof);
    if (it == mDofValues.end())
        return 0;

    return it->second[0].rows();
}

int NuTo::NodeDof::GetDof(Node::eDof rDof, int rComponent) const
{
    const auto& it = mDofNumbers.find(rDof);
    if (it == mDofNumbers.end())
        throw NuTo::MechanicsException(__PRETTY_FUNCTION__, "Cannot access dof type " + Node::DofToString(rDof));

    if (rComponent >= it->second.size())
        throw NuTo::MechanicsException(__PRETTY_FUNCTION__, "Cannot access component " + std::to_string(rComponent) + ". Component " + std::to_string(it->second.size()) + " was requested.");

    return it->second[rComponent];
}

const Eigen::VectorXd& NuTo::NodeDof::Get(Node::eDof rDof, int rTimeDerivative) const
{
    const auto& it = mDofValues.find(rDof);
    if (it == mDofValues.end())
        throw NuTo::MechanicsException(__PRETTY_FUNCTION__, "Cannot access dof type " + Node::DofToString(rDof));

    if (rTimeDerivative >= (int)it->second.size())
        throw NuTo::MechanicsException(__PRETTY_FUNCTION__, "Cannot access time derivative " + std::to_string(rTimeDerivative) + ". This node only has " + std::to_string(it->second.size()) + ".");

    return it->second[rTimeDerivative];
}

void NuTo::NodeDof::Set(Node::eDof rDof, int rTimeDerivative , const Eigen::VectorXd& rValue)
{
    auto it = mDofValues.find(rDof);
    if (it == mDofValues.end())
        throw NuTo::MechanicsException(__PRETTY_FUNCTION__, "Cannot access dof type " + Node::DofToString(rDof));

    if (rTimeDerivative >= (int)it->second.size())
        throw NuTo::MechanicsException(__PRETTY_FUNCTION__, "Cannot access time derivative " + std::to_string(rTimeDerivative) + ". This node only has " + std::to_string(it->second.size()) + ".");

    it->second[rTimeDerivative] = rValue;
}

std::set<NuTo::Node::eDof> NuTo::NodeDof::GetDofTypes() const
{
    std::set<Node::eDof> dofTypes;
    for (auto it : mDofValues)
        dofTypes.insert(it.first);
    return dofTypes;
}

std::string NuTo::NodeDof::GetNodeTypeStr() const
{
    std::stringstream nodeDofype;
    for (auto& it : mDofValues)
    {
        nodeDofype << Node::DofToString(it.first) << ": " << it.second[0].rows() << " dt:" << it.second.size() << "\n";
    }
    return nodeDofype.str();
}

NuTo::NodeBase* NuTo::NodeDof::Clone() const
{
    return new NodeDof(*this);
}

double NuTo::NodeDof::GetDofValueFromVector(
        int rDofNumber,
        const Eigen::VectorXd& rActiveDofValues,
        const Eigen::VectorXd& rDependentDofValues) const
{
    if (rDofNumber < rActiveDofValues.rows())
    {
        return rActiveDofValues(rDofNumber);
    }
    else
    {
        rDofNumber -= rActiveDofValues.rows();
        assert(rDofNumber < rDependentDofValues.rows());
        return rDependentDofValues(rDofNumber);
    }
}

void NuTo::NodeDof::WriteNodeValueToVector(
        int rDofNumber,
        double rDofValue,
        Eigen::VectorXd& rActiveDofValues,
        Eigen::VectorXd& rDependentDofValues) const
{
    if (rDofNumber < rActiveDofValues.rows())
    {
        rActiveDofValues(rDofNumber) = rDofValue;
    }
    else
    {
        rDofNumber -= rActiveDofValues.rows();
        assert(rDofNumber < rDependentDofValues.rows());
        rDependentDofValues(rDofNumber) = rDofValue;
    }
}

#ifdef ENABLE_SERIALIZATION
//! @brief serializes the class
//! @param ar         archive
//! @param version    version
template void NuTo::NodeDof::serialize(boost::archive::binary_oarchive & ar, const unsigned int version);
template void NuTo::NodeDof::serialize(boost::archive::binary_iarchive & ar, const unsigned int version);
template void NuTo::NodeDof::serialize(boost::archive::xml_oarchive & ar, const unsigned int version);
template void NuTo::NodeDof::serialize(boost::archive::xml_iarchive & ar, const unsigned int version);
template void NuTo::NodeDof::serialize(boost::archive::text_oarchive & ar, const unsigned int version);
template void NuTo::NodeDof::serialize(boost::archive::text_iarchive & ar, const unsigned int version);
template<class Archive>
void NuTo::NodeDof::serialize(Archive & ar, const unsigned int version)
{
#ifdef DEBUG_SERIALIZATION
    std::cout << "start serialize NodeDof" << "\n";
#endif
    ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(NodeBase)
       & BOOST_SERIALIZATION_NVP(mDofValues)
       & BOOST_SERIALIZATION_NVP(mDofNumbers);
#ifdef DEBUG_SERIALIZATION
    std::cout << "finish serialize NodeDof \n";
#endif
}
BOOST_CLASS_EXPORT_IMPLEMENT(NuTo::NodeDof)
#endif // ENABLE_SERIALIZATION
