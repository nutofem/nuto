#include "nuto/mechanics/nodes/NodeDof.h"

#ifdef ENABLE_SERIALIZATION
#include "nuto/math/CustomBoostSerializationExtensions.h"
#include <boost/serialization/unordered_map.hpp>
#include <boost/serialization/vector.hpp>
#endif // ENABLE_SERIALIZATION

NuTo::NodeDof::NodeDof(int rNumTimeDerivatives, std::map<Node::eDof, int> rDofDimensions)
     : NodeBase::NodeBase(), mNumTimeDerivatives(rNumTimeDerivatives)
{
    for (auto& it : rDofDimensions)
    {
        Node::eDof dofType = it.first;
        int dofDimension = it.second;

        // allocate global dof numbers
        mDofNumbers[dofType] = Eigen::VectorXi::Zero(dofDimension);


        // allocate dof values
        int numTimeDerivativesForAllocation = mNumTimeDerivatives + 1; // +1 since the 0th time derivative needs to be stored as well

        if (dofType == Node::COORDINATES)
            numTimeDerivativesForAllocation = 1; // no time derivatives for coordinates

        mDofValues[dofType].resize(numTimeDerivativesForAllocation);
        for (int i = 0; i < numTimeDerivativesForAllocation; ++i)
            mDofValues[dofType][i] = Eigen::VectorXd::Zero(dofDimension);
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
        const FullVector<double, Eigen::Dynamic>& rActiveDofValues,
        const FullVector<double, Eigen::Dynamic>& rDependentDofValues)
{
    assert(mNumTimeDerivatives >= rTimeDerivative);
    if (mDofNumbers.find(rDofType) == mDofNumbers.end())
        return;

    for (int i = 0; i < mDofNumbers[rDofType].rows(); ++i)
    {
        double dofValueToSet = GetDofValueFromVector(mDofNumbers[rDofType][i], rActiveDofValues, rDependentDofValues);
        mDofValues[rDofType][rTimeDerivative][i] = dofValueToSet; // this might need some assertions.
    }
}

void NuTo::NodeDof::GetGlobalDofValues(
        int rTimeDerivative,
        Node::eDof rDofType,
        FullVector<double, Eigen::Dynamic>& rActiveDofValues,
        FullVector<double, Eigen::Dynamic>& rDependentDofValues) const
{
    assert(mNumTimeDerivatives >= rTimeDerivative);

    const auto& it = mDofNumbers.find(rDofType);

    if (it == mDofNumbers.end())
        return;

    for (int i = 0; i < it->second.rows(); ++i)
    {
        double dofValueToWrite = mDofValues.at(rDofType)[rTimeDerivative][i];
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



int NuTo::NodeDof::GetNumDofs() const
{
    int numDofs = 0;
    for (auto& it : mDofNumbers)
        numDofs += it.second.rows();
    return numDofs;
}


int NuTo::NodeDof::GetNum(Node::eDof rDof) const
{
    const auto& it = mDofNumbers.find(rDof);
    if (it != mDofNumbers.end())
        return it->second.rows();
    else
        return 0;
}

int NuTo::NodeDof::GetDof(Node::eDof rDof, int rComponent) const
{
    try
    {
        return mDofNumbers.at(rDof)[rComponent];
    }
    catch (...)
    {
        throw NuTo::MechanicsException(__PRETTY_FUNCTION__, "Cannot access dof type " + Node::DofToString(rDof));
    }
}

const Eigen::VectorXd& NuTo::NodeDof::Get(Node::eDof rDof, int rTimeDerivative) const
{
    try
    {
        return mDofValues.at(rDof)[rTimeDerivative];
    }
    catch (...)
    {
        throw NuTo::MechanicsException(__PRETTY_FUNCTION__, "Cannot access dof type " + Node::DofToString(rDof));
    }
}

void NuTo::NodeDof::Set(Node::eDof rDof, int rTimeDerivative , const Eigen::VectorXd& rValue)
{
    if (mDofValues.find(rDof) == mDofValues.end())
        throw NuTo::MechanicsException(__PRETTY_FUNCTION__, "Cannot access dof type " + Node::DofToString(rDof));


    mDofValues[rDof][rTimeDerivative] = rValue;
}

std::set<NuTo::Node::eDof> NuTo::NodeDof::GetDofTypes() const
{
    std::set<Node::eDof> dofTypes;
    for (auto it : mDofNumbers)
        dofTypes.insert(it.first);
    return dofTypes;
}

std::string NuTo::NodeDof::GetNodeTypeStr() const
{
    std::stringstream nodeDofype;

    if (mNumTimeDerivatives > 0)
        nodeDofype << "TimeDerivatives:" << mNumTimeDerivatives << "\n";

    for (auto& it : mDofValues)
        nodeDofype << Node::DofToString(it.first) << ": " << it.second[0].rows() << "\n";

    return nodeDofype.str();

}

NuTo::NodeBase* NuTo::NodeDof::Clone() const
{
    return new NodeDof(*this);
}

double NuTo::NodeDof::GetDofValueFromVector(
        int rDofNumber,
        const FullVector<double, Eigen::Dynamic>& rActiveDofValues,
        const FullVector<double, Eigen::Dynamic>& rDependentDofValues) const
{
    if (rDofNumber < rActiveDofValues.GetNumRows())
    {
        return rActiveDofValues(rDofNumber);
    }
    else
    {
        rDofNumber -= rActiveDofValues.GetNumRows();
        assert(rDofNumber < rDependentDofValues.GetNumRows());
        return rDependentDofValues(rDofNumber);
    }
}

void NuTo::NodeDof::WriteNodeValueToVector(
        int rDofNumber,
        double rDofValue,
        FullVector<double, Eigen::Dynamic>& rActiveDofValues,
        FullVector<double, Eigen::Dynamic>& rDependentDofValues) const
{
    if (rDofNumber < rActiveDofValues.GetNumRows())
    {
        rActiveDofValues(rDofNumber) = rDofValue;
    }
    else
    {
        rDofNumber -= rActiveDofValues.GetNumRows();
        assert(rDofNumber < rDependentDofValues.GetNumRows());
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
       & BOOST_SERIALIZATION_NVP(mNumTimeDerivatives)
       & BOOST_SERIALIZATION_NVP(mDofValues)
       & BOOST_SERIALIZATION_NVP(mDofNumbers);
#ifdef DEBUG_SERIALIZATION
    std::cout << "finish serialize NodeDof \n";
#endif
}
BOOST_CLASS_EXPORT_IMPLEMENT(NuTo::NodeDof)
#endif // ENABLE_SERIALIZATION
