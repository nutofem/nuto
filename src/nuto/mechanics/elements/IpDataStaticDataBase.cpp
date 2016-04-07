// $Id$ 
// IpDataStaticDataBase.cpp
// created Apr 29, 2010 by Joerg F. Unger

#ifdef ENABLE_SERIALIZATION
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/xml_oarchive.hpp>
#include <boost/archive/xml_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <boost/ptr_container/serialize_ptr_vector.hpp>
#endif  // ENABLE_SERIALIZATION

#include "nuto/mechanics/elements/IpDataStaticDataBase.h"
#include "nuto/mechanics/MechanicsException.h"
#include "nuto/mechanics/elements/ElementBase.h"

NuTo::IpDataStaticDataBase::IpDataStaticDataBase() : IpDataBase()
{
}

NuTo::IpDataStaticDataBase::~IpDataStaticDataBase()
{
}

void NuTo::IpDataStaticDataBase::Initialize(const ElementBase* rElement, const ConstitutiveBase* rConstitutive)
{
    mStaticData.clear();

    if (rConstitutive != nullptr)
    {
        auto staticDataRawPtr = rElement->AllocateStaticData(rConstitutive);
        if (staticDataRawPtr != nullptr)
            mStaticData.push_back(staticDataRawPtr);
    }
}

void NuTo::IpDataStaticDataBase::AllocateAdditionalStaticData(int rNumAdditionalStaticData)
{
    if (mStaticData.empty())
        throw MechanicsException(__PRETTY_FUNCTION__, "No static data allocated.");

    for (int i = 0; i < rNumAdditionalStaticData; ++i)
    {
        mStaticData.push_back(GetStaticData(0)->Clone());
    }

}

//! @brief sets the fine scale model (deserialization from a binary file)
void NuTo::IpDataStaticDataBase::SetFineScaleModel(std::string rFileName, double rMacroLength, double rCoordinates[2], std::string rIpName)
{
    if (not mStaticData.empty())
	    mStaticData[0].SetFineScaleModel(rFileName, rMacroLength, rCoordinates, rIpName);
    else
    	throw MechanicsException(__PRETTY_FUNCTION__, "Static data for Ip is not allocated. Either you forgot to assign a material to the ip, or the material law has no static data");
}

//! @brief sets the fine scale parameter for all ips
//! @parameter rName name of the parameter, e.g. YoungsModulus
//! @parameter rParameter value of the parameter
void NuTo::IpDataStaticDataBase::SetFineScaleParameter(const std::string& rName, double rParameter)
{
    if (not mStaticData.empty())
        mStaticData[0].SetFineScaleParameter(rName, rParameter);
    else
    	throw MechanicsException(__PRETTY_FUNCTION__, "Static data for Ip is not allocated. Either you forgot to assign a material to the ip, or the material law has no static data");

}

//! @brief sets the fine scale parameter for all ips
//! @parameter rName name of the parameter, e.g. YoungsModulus
//! @parameter rParameter value of the parameter
void NuTo::IpDataStaticDataBase::SetFineScaleParameter(const std::string& rName, std::string rParameter)
{
    if (not mStaticData.empty())
        mStaticData[0].SetFineScaleParameter(rName, rParameter);
    else
    	throw MechanicsException(__PRETTY_FUNCTION__, "Static data for Ip is not allocated. Either you forgot to assign a material to the ip, or the material law has no static data");

}


//! @brief puts current static data to previous static data, previous to pre-previous, etc.
void NuTo::IpDataStaticDataBase::SaveStaticData()
{
    assert(mStaticData.size() > 1);
    mStaticData.pop_back(); // thanks Sebastian!
    mStaticData.insert(mStaticData.begin(), mStaticData[0].Clone());
}

//! @brief puts current static data to previous static data, previous to pre-previous, etc.
void NuTo::IpDataStaticDataBase::RestoreStaticData()
{
    assert(mStaticData.size() > 1);
    unsigned int numStaticData = mStaticData.size();

    mStaticData.erase(mStaticData.begin());
    // copy the last one
    mStaticData.push_back(mStaticData[numStaticData-2].Clone());
}

int NuTo::IpDataStaticDataBase::GetNumStaticData() const
{
    return mStaticData.size();
}

#ifdef ENABLE_VISUALIZE
//Visualize for all integration points the fine scale structure
void NuTo::IpDataStaticDataBase::VisualizeIpMultiscale(VisualizeUnstructuredGrid& rVisualize,
		const boost::ptr_list<NuTo::VisualizeComponentBase>& rWhat, bool rVisualizeDamage)const
{
    if (not mStaticData.empty())
        mStaticData[0].VisualizeIpMultiscale(rVisualize, rWhat, rVisualizeDamage);
}
#endif


#ifdef ENABLE_SERIALIZATION
// serializes the class
template void NuTo::IpDataStaticDataBase::serialize(boost::archive::binary_oarchive & ar, const unsigned int version);
template void NuTo::IpDataStaticDataBase::serialize(boost::archive::xml_oarchive & ar, const unsigned int version);
template void NuTo::IpDataStaticDataBase::serialize(boost::archive::text_oarchive & ar, const unsigned int version);
template void NuTo::IpDataStaticDataBase::serialize(boost::archive::binary_iarchive & ar, const unsigned int version);
template void NuTo::IpDataStaticDataBase::serialize(boost::archive::xml_iarchive & ar, const unsigned int version);
template void NuTo::IpDataStaticDataBase::serialize(boost::archive::text_iarchive & ar, const unsigned int version);
template<class Archive>
void NuTo::IpDataStaticDataBase::serialize(Archive & ar, const unsigned int version)
{
#ifdef DEBUG_SERIALIZATION
    std::cout << "start serialize IpDataStaticDataBase" << std::endl;
#endif
    ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(IpDataBase)
       & BOOST_SERIALIZATION_NVP(mStaticData);
#ifdef DEBUG_SERIALIZATION
    std::cout << "finish serialize IpDataStaticDataBase" << std::endl;
#endif
}
BOOST_SERIALIZATION_ASSUME_ABSTRACT(NuTo::IpDataStaticDataBase)
#endif // ENABLE_SERIALIZATION
