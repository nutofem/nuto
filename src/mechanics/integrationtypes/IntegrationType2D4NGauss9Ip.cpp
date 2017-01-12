// $Id$

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

#include "mechanics/integrationtypes/IntegrationType2D4NGauss9Ip.h"
#include <assert.h>


//! @brief constructor
NuTo::IntegrationType2D4NGauss9Ip::IntegrationType2D4NGauss9Ip()
{
}

//! @brief returns the local coordinates of an integration point
//! @param rIpNum integration point (counting from zero)
//! @param rCoordinates (result)
void NuTo::IntegrationType2D4NGauss9Ip::GetLocalIntegrationPointCoordinates2D(int rIpNum, double rCoordinates[2])const
{
    assert(rIpNum>=0 && rIpNum<9);
    switch (rIpNum)
    {
    case 0 :
        rCoordinates[0] = -0.774596669241483;
        rCoordinates[1] = -0.774596669241483;
        break;
    case 1 :
        rCoordinates[0] = +0.774596669241483;
        rCoordinates[1] = -0.774596669241483;
        break;
    case 2 :
        rCoordinates[0] = +0.774596669241483;
        rCoordinates[1] = +0.774596669241483;
        break;
    case 3 :
        rCoordinates[0] = -0.774596669241483;
        rCoordinates[1] = +0.774596669241483;
        break;
    case 4 :
        rCoordinates[0] =  0.0;
        rCoordinates[1] = -0.774596669241483;
        break;
    case 5 :
        rCoordinates[0] = +0.774596669241483;
        rCoordinates[1] =  0.0;
        break;
    case 6 :
        rCoordinates[0] =  0.0;
        rCoordinates[1] = +0.774596669241483;
        break;
    case 7 :
        rCoordinates[0] = -0.774596669241483;
        rCoordinates[1] =  0.0;
        break;
    case 8 :
        rCoordinates[0] =  0.0;
        rCoordinates[1] =  0.0;
        break;
    default:
        throw MechanicsException("[NuTo::IntegrationType2D4NGauss9Ip::GetLocalIntegrationPointCoordinates] Ip number out of range.");
    }
}


//! @brief returns the total number of integration points for this integration type
//! @return number of integration points
int NuTo::IntegrationType2D4NGauss9Ip::GetNumIntegrationPoints()const
{
    return 9;
}

//! @brief returns the weight of an integration point
//! @param rIpNum integration point (counting from zero)
//! @return weight of integration points
double NuTo::IntegrationType2D4NGauss9Ip::GetIntegrationPointWeight(int rIpNum)const
{
    assert(rIpNum>=0 && rIpNum<9);
    switch (rIpNum)
    {
    case 0 :
        return 0.308641975; // 5/9 * 5/9
        break;
    case 1 :
        return 0.308641975; // 5/9 * 5/9
        break;
    case 2 :
        return 0.308641975; // 5/9 * 5/9
        break;
    case 3 :
        return 0.308641975; // 5/9 * 5/9
        break;
    case 4 :
        return 0.493827160; // 5/9 * 8/9
        break;
    case 5 :
        return 0.493827160; // 5/9 * 8/9
        break;
    case 6 :
        return 0.493827160; // 5/9 * 8/9
        break;
    case 7 :
        return 0.493827160; // 5/9 * 8/9
        break;
    case 8 :
        return 0.790123456; // 8/9 * 8/9
        break;
    default:
        throw MechanicsException("[NuTo::IntegrationType2D4NGauss9Ip::GetLocalIntegrationPointCoordinates] Ip number out of range.");
    }
}

//! @brief returns a string with the identifier of the integration type
//! @return identifier
std::string NuTo::IntegrationType2D4NGauss9Ip::GetStrIdentifier()const
{
    return GetStrIdentifierStatic();
}

//! @brief returns a string with the identifier of the integration type
//! @return identifier
std::string NuTo::IntegrationType2D4NGauss9Ip::GetStrIdentifierStatic()
{
    return std::string("2D4NGAUSS9IP");
}

#ifdef ENABLE_VISUALIZE
void NuTo::IntegrationType2D4NGauss9Ip::GetVisualizationCells(
    unsigned int& NumVisualizationPoints,
    std::vector<double>& VisualizationPointLocalCoordinates,
    unsigned int& NumVisualizationCells,
    std::vector<NuTo::eCellTypes>& VisualizationCellType,
    std::vector<unsigned int>& VisualizationCellsIncidence,
    std::vector<unsigned int>& VisualizationCellsIP) const
{
    NumVisualizationPoints = 16;

    // first row
    // Point 0
    VisualizationPointLocalCoordinates.push_back(-1);
    VisualizationPointLocalCoordinates.push_back(-1);

    // Point 1
    VisualizationPointLocalCoordinates.push_back(-1./3);
    VisualizationPointLocalCoordinates.push_back(-1);

    // Point 2
    VisualizationPointLocalCoordinates.push_back(+1./3);
    VisualizationPointLocalCoordinates.push_back(-1);

    // Point 3
    VisualizationPointLocalCoordinates.push_back(1);
    VisualizationPointLocalCoordinates.push_back(-1);

    // second row
    // Point 4
    VisualizationPointLocalCoordinates.push_back(-1);
    VisualizationPointLocalCoordinates.push_back(-1./3);

    // Point 5
    VisualizationPointLocalCoordinates.push_back(-1./3);
    VisualizationPointLocalCoordinates.push_back(-1./3);

    // Point 6
    VisualizationPointLocalCoordinates.push_back(+1./3);
    VisualizationPointLocalCoordinates.push_back(-1./3);

    // Point 7
    VisualizationPointLocalCoordinates.push_back(1);
    VisualizationPointLocalCoordinates.push_back(-1./3);

    // third row
    // Point 8
    VisualizationPointLocalCoordinates.push_back(-1);
    VisualizationPointLocalCoordinates.push_back(1./3);

    // Point 9
    VisualizationPointLocalCoordinates.push_back(-1./3);
    VisualizationPointLocalCoordinates.push_back(1./3);

    // Point 10
    VisualizationPointLocalCoordinates.push_back(+1./3);
    VisualizationPointLocalCoordinates.push_back(1./3);

    // Point 11
    VisualizationPointLocalCoordinates.push_back(1);
    VisualizationPointLocalCoordinates.push_back(1./3);

    // fourth row
    // Point 12
    VisualizationPointLocalCoordinates.push_back(-1);
    VisualizationPointLocalCoordinates.push_back(1);

    // Point 13
    VisualizationPointLocalCoordinates.push_back(-1./3);
    VisualizationPointLocalCoordinates.push_back(1);

    // Point 14
    VisualizationPointLocalCoordinates.push_back(+1./3);
    VisualizationPointLocalCoordinates.push_back(1);

    // Point 15
    VisualizationPointLocalCoordinates.push_back(1);
    VisualizationPointLocalCoordinates.push_back(1);

    NumVisualizationCells = 9;

    // cell 0
    VisualizationCellType.push_back(NuTo::eCellTypes::QUAD);
    VisualizationCellsIncidence.push_back(0);
    VisualizationCellsIncidence.push_back(1);
    VisualizationCellsIncidence.push_back(5);
    VisualizationCellsIncidence.push_back(4);
    VisualizationCellsIP.push_back(0);

    // cell 4
    VisualizationCellType.push_back(NuTo::eCellTypes::QUAD);
    VisualizationCellsIncidence.push_back(1);
    VisualizationCellsIncidence.push_back(2);
    VisualizationCellsIncidence.push_back(6);
    VisualizationCellsIncidence.push_back(5);
    VisualizationCellsIP.push_back(4);

    // cell 1
    VisualizationCellType.push_back(NuTo::eCellTypes::QUAD);
    VisualizationCellsIncidence.push_back(2);
    VisualizationCellsIncidence.push_back(3);
    VisualizationCellsIncidence.push_back(7);
    VisualizationCellsIncidence.push_back(6);
    VisualizationCellsIP.push_back(1);

    // cell 7
    VisualizationCellType.push_back(NuTo::eCellTypes::QUAD);
    VisualizationCellsIncidence.push_back(4);
    VisualizationCellsIncidence.push_back(5);
    VisualizationCellsIncidence.push_back(9);
    VisualizationCellsIncidence.push_back(8);
    VisualizationCellsIP.push_back(7);

    // cell 8
    VisualizationCellType.push_back(NuTo::eCellTypes::QUAD);
    VisualizationCellsIncidence.push_back(5);
    VisualizationCellsIncidence.push_back(6);
    VisualizationCellsIncidence.push_back(10);
    VisualizationCellsIncidence.push_back(9);
    VisualizationCellsIP.push_back(8);

    // cell 5
    VisualizationCellType.push_back(NuTo::eCellTypes::QUAD);
    VisualizationCellsIncidence.push_back(6);
    VisualizationCellsIncidence.push_back(7);
    VisualizationCellsIncidence.push_back(11);
    VisualizationCellsIncidence.push_back(10);
    VisualizationCellsIP.push_back(5);

    // cell 3
    VisualizationCellType.push_back(NuTo::eCellTypes::QUAD);
    VisualizationCellsIncidence.push_back(8);
    VisualizationCellsIncidence.push_back(9);
    VisualizationCellsIncidence.push_back(13);
    VisualizationCellsIncidence.push_back(12);
    VisualizationCellsIP.push_back(3);

    // cell 6
    VisualizationCellType.push_back(NuTo::eCellTypes::QUAD);
    VisualizationCellsIncidence.push_back(9);
    VisualizationCellsIncidence.push_back(10);
    VisualizationCellsIncidence.push_back(14);
    VisualizationCellsIncidence.push_back(13);
    VisualizationCellsIP.push_back(6);

    // cell 2
    VisualizationCellType.push_back(NuTo::eCellTypes::QUAD);
    VisualizationCellsIncidence.push_back(10);
    VisualizationCellsIncidence.push_back(11);
    VisualizationCellsIncidence.push_back(15);
    VisualizationCellsIncidence.push_back(14);
    VisualizationCellsIP.push_back(2);

}
#endif // ENABLE_VISUALIZE

#ifdef ENABLE_SERIALIZATION
// serializes the class
template void NuTo::IntegrationType2D4NGauss9Ip::serialize(boost::archive::binary_oarchive & ar, const unsigned int version);
template void NuTo::IntegrationType2D4NGauss9Ip::serialize(boost::archive::xml_oarchive & ar, const unsigned int version);
template void NuTo::IntegrationType2D4NGauss9Ip::serialize(boost::archive::text_oarchive & ar, const unsigned int version);
template void NuTo::IntegrationType2D4NGauss9Ip::serialize(boost::archive::binary_iarchive & ar, const unsigned int version);
template void NuTo::IntegrationType2D4NGauss9Ip::serialize(boost::archive::xml_iarchive & ar, const unsigned int version);
template void NuTo::IntegrationType2D4NGauss9Ip::serialize(boost::archive::text_iarchive & ar, const unsigned int version);
template<class Archive>
void NuTo::IntegrationType2D4NGauss9Ip::serialize(Archive & ar, const unsigned int version)
{
#ifdef DEBUG_SERIALIZATION
    std::cout << "start serialize IntegrationType2D4NGauss9Ip" << std::endl;
#endif
    ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(IntegrationType2D);
#ifdef DEBUG_SERIALIZATION
    std::cout << "finish serialize IntegrationType2D4NGauss9Ip" << std::endl;
#endif
}
BOOST_CLASS_EXPORT_IMPLEMENT(NuTo::IntegrationType2D4NGauss9Ip)
#endif // ENABLE_SERIALIZATION
