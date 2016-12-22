// $Id: IntegrationType2D3NGauss13Ip.cpp 309 2010-09-22 22:21:24Z unger3 $

#ifdef ENABLE_SERIALIZATION
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/xml_oarchive.hpp>
#include <boost/archive/xml_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#endif //ENABLE_SERIALIZATION

#ifdef ENABLE_VISUALIZE
#include "visualize/VisualizeEnum.h"
#endif // ENABLE_VISUALIZE

#include "mechanics/integrationtypes/IntegrationType2D3NGauss13Ip.h"
#include <assert.h>


//! @brief constructor
NuTo::IntegrationType2D3NGauss13Ip::IntegrationType2D3NGauss13Ip()
{
}

//! @brief returns the local coordinates of an integration point
//! @param rIpNum integration point (counting from zero)
//! @param rCoordinates (result)
void NuTo::IntegrationType2D3NGauss13Ip::GetLocalIntegrationPointCoordinates2D(int rIpNum, double rCoordinates[2])const
{
    assert(rIpNum>=0 && rIpNum<13);
    switch (rIpNum)
    {
    case  0  :
        rCoordinates[0]=  0.333333333333 ;
        rCoordinates[1]=  0.333333333333 ;
        break;
    case  1  :
        rCoordinates[0]=  0.0523383720927 ;
        rCoordinates[1]=  0.473830813954 ;
        break;
    case  2  :
        rCoordinates[0]=  0.473830813954 ;
        rCoordinates[1]=  0.473830813954 ;
        break;
    case  3  :
        rCoordinates[0]=  0.473830813954 ;
        rCoordinates[1]=  0.0523383720927 ;
        break;
    case  4  :
        rCoordinates[0]=  0.655764660738 ;
        rCoordinates[1]=  0.172117669631 ;
        break;
    case  5  :
        rCoordinates[0]=  0.172117669631 ;
        rCoordinates[1]=  0.172117669631 ;
        break;
    case  6  :
        rCoordinates[0]=  0.172117669631 ;
        rCoordinates[1]=  0.655764660738 ;
        break;
    case  7  :
        rCoordinates[0]=  0.0 ;
        rCoordinates[1]=  0.865307354083 ;
        break;
    case  8  :
        rCoordinates[0]=  0.0 ;
        rCoordinates[1]=  0.134692645917 ;
        break;
    case  9  :
        rCoordinates[0]=  0.865307354083 ;
        rCoordinates[1]=  0.0 ;
        break;
    case  10  :
        rCoordinates[0]=  0.865307354083 ;
        rCoordinates[1]=  0.134692645917 ;
        break;
    case  11  :
        rCoordinates[0]=  0.134692645917 ;
        rCoordinates[1]=  0.0 ;
        break;
    case  12  :
        rCoordinates[0]=  0.134692645917 ;
        rCoordinates[1]=  0.865307354083 ;
        break;
    default:
        throw MechanicsException("[NuTo::IntegrationType2D3NGauss13Ip::GetLocalIntegrationPointCoordinates] Ip number out of range.");
    }
}


//! @brief returns the total number of integration points for this integration type
//! @return number of integration points
int NuTo::IntegrationType2D3NGauss13Ip::GetNumIntegrationPoints()const
{
    return 13;
}

//! @brief returns the weight of an integration point
//! @param rIpNum integration point (counting from zero)
//! @return weight of integration points
double NuTo::IntegrationType2D3NGauss13Ip::GetIntegrationPointWeight(int rIpNum)const
{
    assert(rIpNum>=0 && rIpNum<13);
    switch (rIpNum)
    {
    case  0  :
        return  0.0763544833942 ;
        break;
    case  1  :
        return  0.0490679340394 ;
        break;
    case  2  :
        return  0.0490679340394 ;
        break;
    case  3  :
        return  0.0490679340394 ;
        break;
    case  4  :
        return  0.0647842146403 ;
        break;
    case  5  :
        return  0.0647842146403 ;
        break;
    case  6  :
        return  0.0647842146403 ;
        break;
    case  7  :
        return  0.0136815117611 ;
        break;
    case  8  :
        return  0.0136815117611 ;
        break;
    case  9  :
        return  0.0136815117611 ;
        break;
    case  10  :
        return  0.0136815117611 ;
        break;
    case  11  :
        return  0.0136815117611 ;
        break;
    case  12  :
        return  0.0136815117611 ;
        break;
    default:
        throw MechanicsException("[NuTo::IntegrationType2D3NGauss13Ip::GetLocalIntegrationPointCoordinates] Ip number out of range.");
    }
}

//! @brief returns a string with the identifier of the integration type
//! @return identifier
std::string NuTo::IntegrationType2D3NGauss13Ip::GetStrIdentifier()const
{
    return GetStrIdentifierStatic();
}

//! @brief returns a string with the identifier of the integration type
//! @return identifier
std::string NuTo::IntegrationType2D3NGauss13Ip::GetStrIdentifierStatic()
{
    return std::string("2D3NGAUSS13IP");
}

#ifdef ENABLE_VISUALIZE
void NuTo::IntegrationType2D3NGauss13Ip::GetVisualizationCells(
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
    VisualizationCellsIP.push_back(4);

    // cell 2
    VisualizationCellType.push_back(NuTo::eCellTypes::QUAD);
    VisualizationCellsIncidence.push_back(4);
    VisualizationCellsIncidence.push_back(5);
    VisualizationCellsIncidence.push_back(6);
    VisualizationCellsIncidence.push_back(3);
    VisualizationCellsIP.push_back(6);
}
#endif // ENABLE_VISUALIZE

#ifdef ENABLE_SERIALIZATION
// serializes the class
template void NuTo::IntegrationType2D3NGauss13Ip::serialize(boost::archive::binary_oarchive & ar, const unsigned int version);
template void NuTo::IntegrationType2D3NGauss13Ip::serialize(boost::archive::xml_oarchive & ar, const unsigned int version);
template void NuTo::IntegrationType2D3NGauss13Ip::serialize(boost::archive::text_oarchive & ar, const unsigned int version);
template void NuTo::IntegrationType2D3NGauss13Ip::serialize(boost::archive::binary_iarchive & ar, const unsigned int version);
template void NuTo::IntegrationType2D3NGauss13Ip::serialize(boost::archive::xml_iarchive & ar, const unsigned int version);
template void NuTo::IntegrationType2D3NGauss13Ip::serialize(boost::archive::text_iarchive & ar, const unsigned int version);
template<class Archive>
void NuTo::IntegrationType2D3NGauss13Ip::serialize(Archive & ar, const unsigned int version)
{
#ifdef DEBUG_SERIALIZATION
    std::cout << "start serialize IntegrationType2D3NGauss13Ip" << std::endl;
#endif
    ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(IntegrationType2D);
#ifdef DEBUG_SERIALIZATION
    std::cout << "finish serialize IntegrationType2D3NGauss13Ip" << std::endl;
#endif
}
BOOST_CLASS_EXPORT_IMPLEMENT(NuTo::IntegrationType2D3NGauss13Ip)
#endif // ENABLE_SERIALIZATION
