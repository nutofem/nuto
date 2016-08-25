// $Id: IntegrationType2D3NGauss16Ip.cpp 309 2010-09-22 22:21:24Z unger3 $

#ifdef ENABLE_SERIALIZATION
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/xml_oarchive.hpp>
#include <boost/archive/xml_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#endif //ENABLE_SERIALIZATION

#ifdef ENABLE_VISUALIZE
#include "nuto/visualize/VisualizeEnum.h"
#endif // ENABLE_VISUALIZE

#include "nuto/mechanics/integrationtypes/IntegrationType2D3NGauss16Ip.h"
#include <assert.h>


//! @brief constructor
NuTo::IntegrationType2D3NGauss16Ip::IntegrationType2D3NGauss16Ip()
{
}

//! @brief returns the local coordinates of an integration point
//! @param rIpNum integration point (counting from zero)
//! @param rCoordinates (result)
void NuTo::IntegrationType2D3NGauss16Ip::GetLocalIntegrationPointCoordinates2D(int rIpNum, double rCoordinates[2])const
{
    assert(rIpNum>=0 && rIpNum<16);
    switch (rIpNum)
    {
    case  0  :
        rCoordinates[0]=  0.333333333333 ;
        rCoordinates[1]=  0.333333333333 ;
        break;
    case  1  :
        rCoordinates[0]=  0.0814148234146 ;
        rCoordinates[1]=  0.459292588293 ;
        break;
    case  2  :
        rCoordinates[0]=  0.459292588293 ;
        rCoordinates[1]=  0.459292588293 ;
        break;
    case  3  :
        rCoordinates[0]=  0.459292588293 ;
        rCoordinates[1]=  0.0814148234146 ;
        break;
    case  4  :
        rCoordinates[0]=  0.898905543366 ;
        rCoordinates[1]=  0.050547228317 ;
        break;
    case  5  :
        rCoordinates[0]=  0.050547228317 ;
        rCoordinates[1]=  0.050547228317 ;
        break;
    case  6  :
        rCoordinates[0]=  0.050547228317 ;
        rCoordinates[1]=  0.898905543366 ;
        break;
    case  7  :
        rCoordinates[0]=  0.658861384496 ;
        rCoordinates[1]=  0.170569307752 ;
        break;
    case  8  :
        rCoordinates[0]=  0.170569307752 ;
        rCoordinates[1]=  0.170569307752 ;
        break;
    case  9  :
        rCoordinates[0]=  0.170569307752 ;
        rCoordinates[1]=  0.658861384496 ;
        break;
    case  10  :
        rCoordinates[0]=  0.00839477740996 ;
        rCoordinates[1]=  0.728492392955 ;
        break;
    case  11  :
        rCoordinates[0]=  0.00839477740996 ;
        rCoordinates[1]=  0.263112829635 ;
        break;
    case  12  :
        rCoordinates[0]=  0.728492392955 ;
        rCoordinates[1]=  0.00839477740996 ;
        break;
    case  13  :
        rCoordinates[0]=  0.728492392955 ;
        rCoordinates[1]=  0.263112829635 ;
        break;
    case  14  :
        rCoordinates[0]=  0.263112829635 ;
        rCoordinates[1]=  0.00839477740996 ;
        break;
    case  15  :
        rCoordinates[0]=  0.263112829635 ;
        rCoordinates[1]=  0.728492392955 ;
        break;
    default:
        throw MechanicsException("[NuTo::IntegrationType2D3NGauss16Ip::GetLocalIntegrationPointCoordinates] Ip number out of range.");
    }
}


//! @brief returns the total number of integration points for this integration type
//! @return number of integration points
int NuTo::IntegrationType2D3NGauss16Ip::GetNumIntegrationPoints()const
{
    return 16;
}

//! @brief returns the weight of an integration point
//! @param rIpNum integration point (counting from zero)
//! @return weight of integration points
double NuTo::IntegrationType2D3NGauss16Ip::GetIntegrationPointWeight(int rIpNum)const
{
    assert(rIpNum>=0 && rIpNum<16);
    switch (rIpNum)
    {
    case  0  :
        return  0.0721578038389 ;
        break;
    case  1  :
        return  0.0475458171336 ;
        break;
    case  2  :
        return  0.0475458171336 ;
        break;
    case  3  :
        return  0.0475458171336 ;
        break;
    case  4  :
        return  0.0162292488116 ;
        break;
    case  5  :
        return  0.0162292488116 ;
        break;
    case  6  :
        return  0.0162292488116 ;
        break;
    case  7  :
        return  0.0516086852674 ;
        break;
    case  8  :
        return  0.0516086852674 ;
        break;
    case  9  :
        return  0.0516086852674 ;
        break;
    case  10  :
        return  0.0136151570872 ;
        break;
    case  11  :
        return  0.0136151570872 ;
        break;
    case  12  :
        return  0.0136151570872 ;
        break;
    case  13  :
        return  0.0136151570872 ;
        break;
    case  14  :
        return  0.0136151570872 ;
        break;
    case  15  :
        return  0.0136151570872 ;
        break;
    default:
        throw MechanicsException("[NuTo::IntegrationType2D3NGauss16Ip::GetLocalIntegrationPointCoordinates] Ip number out of range.");
    }
}

//! @brief returns a string with the identifier of the integration type
//! @return identifier
std::string NuTo::IntegrationType2D3NGauss16Ip::GetStrIdentifier()const
{
    return GetStrIdentifierStatic();
}

//! @brief returns a string with the identifier of the integration type
//! @return identifier
std::string NuTo::IntegrationType2D3NGauss16Ip::GetStrIdentifierStatic()
{
    return std::string("2D3NGAUSS16IP");
}

#ifdef ENABLE_VISUALIZE
void NuTo::IntegrationType2D3NGauss16Ip::GetVisualizationCells(
    unsigned int& NumVisualizationPoints,
    std::vector<double>& VisualizationPointLocalCoordinates,
    unsigned int& NumVisualizationCells,
    std::vector<NuTo::eCellTypes>& VisualizationCellType,
    std::vector<unsigned int>& VisualizationCellsIncidence,
    std::vector<unsigned int>& VisualizationCellsIP) const
{
    NumVisualizationPoints = 7;

    // Point 0
    VisualizationPointLocalCoordinates.push_back(0);
    VisualizationPointLocalCoordinates.push_back(0);

    // Point 1
    VisualizationPointLocalCoordinates.push_back(0.5);
    VisualizationPointLocalCoordinates.push_back(0);

    // Point 2
    VisualizationPointLocalCoordinates.push_back(1);
    VisualizationPointLocalCoordinates.push_back(0);

    // Point 3
    VisualizationPointLocalCoordinates.push_back(0);
    VisualizationPointLocalCoordinates.push_back(0.5);

    // Point 4
    VisualizationPointLocalCoordinates.push_back(1./3.);
    VisualizationPointLocalCoordinates.push_back(1./3.);

    // Point 5
    VisualizationPointLocalCoordinates.push_back(0.5);
    VisualizationPointLocalCoordinates.push_back(0.5);

    // Point 6
    VisualizationPointLocalCoordinates.push_back(0);
    VisualizationPointLocalCoordinates.push_back(1);

    NumVisualizationCells = 3;

    // cell 0
    VisualizationCellType.push_back(NuTo::eCellTypes::QUAD);
    VisualizationCellsIncidence.push_back(0);
    VisualizationCellsIncidence.push_back(1);
    VisualizationCellsIncidence.push_back(4);
    VisualizationCellsIncidence.push_back(3);
    VisualizationCellsIP.push_back(5);

    // cell 1
    VisualizationCellType.push_back(NuTo::eCellTypes::QUAD);
    VisualizationCellsIncidence.push_back(1);
    VisualizationCellsIncidence.push_back(2);
    VisualizationCellsIncidence.push_back(5);
    VisualizationCellsIncidence.push_back(4);
    VisualizationCellsIP.push_back(6);

    // cell 2
    VisualizationCellType.push_back(NuTo::eCellTypes::QUAD);
    VisualizationCellsIncidence.push_back(4);
    VisualizationCellsIncidence.push_back(5);
    VisualizationCellsIncidence.push_back(6);
    VisualizationCellsIncidence.push_back(3);
    VisualizationCellsIP.push_back(4);
}
#endif // ENABLE_VISUALIZE

#ifdef ENABLE_SERIALIZATION
// serializes the class
template void NuTo::IntegrationType2D3NGauss16Ip::serialize(boost::archive::binary_oarchive & ar, const unsigned int version);
template void NuTo::IntegrationType2D3NGauss16Ip::serialize(boost::archive::xml_oarchive & ar, const unsigned int version);
template void NuTo::IntegrationType2D3NGauss16Ip::serialize(boost::archive::text_oarchive & ar, const unsigned int version);
template void NuTo::IntegrationType2D3NGauss16Ip::serialize(boost::archive::binary_iarchive & ar, const unsigned int version);
template void NuTo::IntegrationType2D3NGauss16Ip::serialize(boost::archive::xml_iarchive & ar, const unsigned int version);
template void NuTo::IntegrationType2D3NGauss16Ip::serialize(boost::archive::text_iarchive & ar, const unsigned int version);
template<class Archive>
void NuTo::IntegrationType2D3NGauss16Ip::serialize(Archive & ar, const unsigned int version)
{
#ifdef DEBUG_SERIALIZATION
    std::cout << "start serialize IntegrationType2D3NGauss16Ip" << std::endl;
#endif
    ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(IntegrationType2D);
#ifdef DEBUG_SERIALIZATION
    std::cout << "finish serialize IntegrationType2D3NGauss16Ip" << std::endl;
#endif
}
BOOST_CLASS_EXPORT_IMPLEMENT(NuTo::IntegrationType2D3NGauss16Ip)
#endif // ENABLE_SERIALIZATION
