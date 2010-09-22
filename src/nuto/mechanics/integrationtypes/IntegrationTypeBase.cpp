#ifdef ENABLE_SERIALIZATION
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/xml_oarchive.hpp>
#include <boost/archive/xml_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <boost/serialization/vector.hpp>
#endif // ENABLE_SERIALIZATION

#include "nuto/mechanics/integrationtypes/IntegrationTypeBase.h"
//! @brief constructor
NuTo::IntegrationTypeBase::IntegrationTypeBase()
{}

//! @brief info about the integration type
//! @param rVerboseLevel determines how detailed the information is
void NuTo::IntegrationTypeBase::Info(int rVerboseLevel)const
{
    std::cout << GetStrIdentifier() << std::endl;
    if (rVerboseLevel>2)
    {
        double localCoord[3];
        for (int count=0; count<GetNumIntegrationPoints(); count++)
        {
            std::cout << "    IP " << count << " weight " << GetIntegrationPointWeight(count) << std::endl;
            std::cout << "        coordinates " ;
            switch (GetCoordinateDimension())
            {
            case 1:
                GetLocalIntegrationPointCoordinates1D(count,localCoord[0]);
                std::cout << localCoord[0] << std::endl;
                break;
            case 2:
                GetLocalIntegrationPointCoordinates2D(count,localCoord);
                std::cout << localCoord[0] << localCoord[1] << std::endl;
                break;
            case 3:
                GetLocalIntegrationPointCoordinates3D(count,localCoord);
                std::cout << localCoord[0] << localCoord[1] << localCoord[2] << std::endl;
                break;
            default:
                throw MechanicsException("[NuTo::IntegrationTypeBase::Info] Invalid dimension of integration point coordinates.");
            }
        }
    }
}

//! @brief returns the local coordinates of an integration point
//! @param rIpNum integration point (counting from zero)
//! @param rCoordinates (result)
void NuTo::IntegrationTypeBase::GetLocalIntegrationPointCoordinates1D(int rIpNum, double& rCoordinates)const
{
    throw MechanicsException("[NuTo::IntegrationTypeBase::GetLocalIntegrationPointCoordinates] This integration type does not support 1D coordinates.");
}

//! @brief returns the local coordinates of an integration point
//! @param rIpNum integration point (counting from zero)
//! @param rCoordinates (result)
void NuTo::IntegrationTypeBase::GetLocalIntegrationPointCoordinates2D(int rIpNum, double rCoordinates[2])const
{
    throw MechanicsException("[NuTo::IntegrationTypeBase::GetLocalIntegrationPointCoordinates] This integration type does not support 2D coordinates.");
}

//! @brief returns the local coordinates of an integration point
//! @param rIpNum integration point (counting from zero)
//! @param rCoordinates (result)
void NuTo::IntegrationTypeBase::GetLocalIntegrationPointCoordinates3D(int rIpNum, double rCoordinates[3])const
{
    throw MechanicsException("[NuTo::IntegrationTypeBase::GetLocalIntegrationPointCoordinates] This integration type does not support 3D coordinates.");
}

#ifdef ENABLE_SERIALIZATION
// serializes the class
template void NuTo::IntegrationTypeBase::serialize(boost::archive::binary_oarchive & ar, const unsigned int version);
template void NuTo::IntegrationTypeBase::serialize(boost::archive::xml_oarchive & ar, const unsigned int version);
template void NuTo::IntegrationTypeBase::serialize(boost::archive::text_oarchive & ar, const unsigned int version);
template void NuTo::IntegrationTypeBase::serialize(boost::archive::binary_iarchive & ar, const unsigned int version);
template void NuTo::IntegrationTypeBase::serialize(boost::archive::xml_iarchive & ar, const unsigned int version);
template void NuTo::IntegrationTypeBase::serialize(boost::archive::text_iarchive & ar, const unsigned int version);
template<class Archive>
void NuTo::IntegrationTypeBase::serialize(Archive & ar, const unsigned int version)
{
#ifdef DEBUG_SERIALIZATION
    std::cout << "start serialize IntegrationTypeBase" << std::endl;
#endif
#ifdef DEBUG_SERIALIZATION
    std::cout << "finish serialize IntegrationTypeBase" << std::endl;
#endif
}
BOOST_CLASS_EXPORT_IMPLEMENT(NuTo::IntegrationTypeBase)
BOOST_SERIALIZATION_ASSUME_ABSTRACT(NuTo::IntegrationTypeBase)
#endif // ENABLE_SERIALIZATION
