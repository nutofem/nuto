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

#include "mechanics/integrationtypes/IntegrationType2D3NGauss12Ip.h"
#include <assert.h>

//! @brief constructor
NuTo::IntegrationType2D3NGauss12Ip::IntegrationType2D3NGauss12Ip()
{
}

//! @brief returns the local coordinates of an integration point
//! @param rIpNum integration point (counting from zero)
//! @param rCoordinates (result)
void NuTo::IntegrationType2D3NGauss12Ip::GetLocalIntegrationPointCoordinates2D(int rIpNum, double rCoordinates[2]) const
{
    assert(rIpNum >= 0 && rIpNum < 12);

    const double a = 0.063089104491502;
    const double b = 0.249286745170910;
    const double c = 0.310352451033785;
    const double d = 0.053145049844816;

    switch (rIpNum)
    {
    case 0:
        rCoordinates[0] = a;
        rCoordinates[1] = a;
        break;
    case 1:
        rCoordinates[0] = 1 - 2 * a;
        rCoordinates[1] = a;
        break;
    case 2:
        rCoordinates[0] = a;
        rCoordinates[1] = 1 - 2 * a;
        break;
    case 3:
        rCoordinates[0] = b;
        rCoordinates[1] = b;
        break;
    case 4:
        rCoordinates[0] = 1 - 2 * b;
        rCoordinates[1] = b;
        break;
    case 5:
        rCoordinates[0] = b;
        rCoordinates[1] = 1 - 2 * b;
        break;
    case 6:
        rCoordinates[0] = c;
        rCoordinates[1] = d;
        break;
    case 7:
        rCoordinates[0] = d;
        rCoordinates[1] = c;
        break;
    case 8:
        rCoordinates[0] = 1 - c - d;
        rCoordinates[1] = c;
        break;
    case 9:
        rCoordinates[0] = 1 - c - d;
        rCoordinates[1] = d;
        break;
    case 10:
        rCoordinates[0] = c;
        rCoordinates[1] = 1 - c - d;
        break;
    case 11:
        rCoordinates[0] = d;
        rCoordinates[1] = 1 - c - d;
        break;
    default:
        throw MechanicsException("[NuTo::IntegrationType2D3NGauss12Ip::GetLocalIntegrationPointCoordinates] Ip number out of range.");
    }
}

//! @brief returns the total number of integration points for this integration type
//! @return number of integration points
int NuTo::IntegrationType2D3NGauss12Ip::GetNumIntegrationPoints() const
{
    return 12;
}

//! @brief returns the weight of an integration point
//! @param rIpNum integration point (counting from zero)
//! @return weight of integration points
double NuTo::IntegrationType2D3NGauss12Ip::GetIntegrationPointWeight(int rIpNum) const
{
    assert(rIpNum >= 0 && rIpNum < 12);
    const double e = 0.025422453185103;
    const double f = 0.058393137863189;
    const double g = 0.041425537809187;
    switch (rIpNum)
    {
    case 0:
        return e;
        break;
    case 1:
        return e;
        break;
    case 2:
        return e;
        break;
    case 3:
        return f;
        break;
    case 4:
        return f;
        break;
    case 5:
        return f;
        break;
    case 6:
        return g;
        break;
    case 7:
        return g;
        break;
    case 8:
        return g;
        break;
    case 9:
        return g;
        break;
    case 10:
        return g;
        break;
    case 11:
        return g;
        break;
    default:
        throw MechanicsException("[NuTo::IntegrationType2D3NGauss12Ip::GetLocalIntegrationPointCoordinates] Ip number out of range.");
    }
}

//! @brief returns a string with the identifier of the integration type
//! @return identifier
std::string NuTo::IntegrationType2D3NGauss12Ip::GetStrIdentifier() const
{
    return GetStrIdentifierStatic();
}

//! @brief returns a string with the identifier of the integration type
//! @return identifier
std::string NuTo::IntegrationType2D3NGauss12Ip::GetStrIdentifierStatic()
{
    return std::string("2D3NGAUSS12IP");
}

#ifdef ENABLE_VISUALIZE
void NuTo::IntegrationType2D3NGauss12Ip::GetVisualizationCells(unsigned int& NumVisualizationPoints, std::vector<double>& VisualizationPointLocalCoordinates, unsigned int& NumVisualizationCells, std::vector<NuTo::eCellTypes>& VisualizationCellType,
        std::vector<unsigned int>& VisualizationCellsIncidence, std::vector<unsigned int>& VisualizationCellsIP) const
{

    // only 3 integration points (3,4,5) are visualised. TODO: Voronoi decomposition + triangulation for proper visualisation

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
    VisualizationPointLocalCoordinates.push_back(1. / 3.);
    VisualizationPointLocalCoordinates.push_back(1. / 3.);

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
    VisualizationCellsIP.push_back(3);

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
    VisualizationCellsIP.push_back(5);
}
#endif // ENABLE_VISUALIZE

#ifdef ENABLE_SERIALIZATION
// serializes the class
template void NuTo::IntegrationType2D3NGauss12Ip::serialize(boost::archive::binary_oarchive & ar, const unsigned int version);
template void NuTo::IntegrationType2D3NGauss12Ip::serialize(boost::archive::xml_oarchive & ar, const unsigned int version);
template void NuTo::IntegrationType2D3NGauss12Ip::serialize(boost::archive::text_oarchive & ar, const unsigned int version);
template void NuTo::IntegrationType2D3NGauss12Ip::serialize(boost::archive::binary_iarchive & ar, const unsigned int version);
template void NuTo::IntegrationType2D3NGauss12Ip::serialize(boost::archive::xml_iarchive & ar, const unsigned int version);
template void NuTo::IntegrationType2D3NGauss12Ip::serialize(boost::archive::text_iarchive & ar, const unsigned int version);
template<class Archive>
void NuTo::IntegrationType2D3NGauss12Ip::serialize(Archive & ar, const unsigned int version)
{
#ifdef DEBUG_SERIALIZATION
    std::cout << "start serialize IntegrationType2D3NGauss12Ip" << std::endl;
#endif
    ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(IntegrationType2D);
#ifdef DEBUG_SERIALIZATION
    std::cout << "finish serialize IntegrationType2D3NGauss12Ip" << std::endl;
#endif
}
BOOST_CLASS_EXPORT_IMPLEMENT(NuTo::IntegrationType2D3NGauss12Ip)
#endif // ENABLE_SERIALIZATION
