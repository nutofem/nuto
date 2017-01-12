// $Id$

#ifdef ENABLE_SERIALIZATION
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/xml_oarchive.hpp>
#include <boost/archive/xml_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <boost/serialization/vector.hpp>
#endif // ENABLE_SERIALIZATION

#include "mechanics/integrationtypes/IntegrationTypeBase.h"

NuTo::IntegrationTypeBase::IntegrationTypeBase()
{}

void NuTo::IntegrationTypeBase::Info(int rVerboseLevel)const
{
    std::cout << IntegrationTypeToString(GetEnumType()) << std::endl;
    if (rVerboseLevel>2)
    {
        for (int count=0; count<GetNumIntegrationPoints(); count++)
        {
            std::cout << "    IP " << count << " weight " << GetIntegrationPointWeight(count) << std::endl;
            std::cout << "        coordinates " ;
            std::cout << GetLocalIntegrationPointCoordinates(count);
        }
    }
}

void NuTo::IntegrationTypeBase::AddIntegrationPoints(std::vector< std::vector<double> > & rArea, const unsigned short rOrder)
{
    throw MechanicsException(__PRETTY_FUNCTION__, "Cannot add an IP to integration type " + IntegrationTypeToString(GetEnumType()) + ".");
}

void NuTo::IntegrationTypeBase::AddIntegrationPoint(const IntegrationPointBase & rIp)
{
    throw MechanicsException(__PRETTY_FUNCTION__, "Cannot add an IP to integration type " + IntegrationTypeToString(GetEnumType()) + ".");
}

void NuTo::IntegrationTypeBase::DeleteIntegrationPoint(const int rIpNum)
{
    throw MechanicsException(__PRETTY_FUNCTION__, "Cannot delete an IP to integration type " + IntegrationTypeToString(GetEnumType()) + ".");
}

#ifdef ENABLE_SERIALIZATION
BOOST_CLASS_EXPORT_IMPLEMENT(NuTo::IntegrationTypeBase)
BOOST_SERIALIZATION_ASSUME_ABSTRACT(NuTo::IntegrationTypeBase)
#endif // ENABLE_SERIALIZATION
