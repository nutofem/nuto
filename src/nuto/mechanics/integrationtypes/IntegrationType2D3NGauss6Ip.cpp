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

#include "nuto/mechanics/integrationtypes/IntegrationType2D3NGauss6Ip.h"
#include <assert.h>

NuTo::IntegrationType2D3NGauss6Ip::IntegrationType2D3NGauss6Ip() {}

void NuTo::IntegrationType2D3NGauss6Ip::GetLocalIntegrationPointCoordinates2D(int rIpNum, double rCoordinates[2]) const
{
    assert(rIpNum >= 0 && rIpNum < 6);

    const double a = 0.445948490915965;
    const double b = 0.091576213509771;

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
    default:
        throw MechanicsException(__PRETTY_FUNCTION__, "Ip number out of range.");
    }
}


unsigned int NuTo::IntegrationType2D3NGauss6Ip::GetNumIntegrationPoints() const
{
    return 6;
}


double NuTo::IntegrationType2D3NGauss6Ip::GetIntegrationPointWeight(int rIpNum) const
{
    const double c = 0.111690794839005;
    const double d = 0.054975871827661;

    assert(rIpNum >= 0 && rIpNum < 6);
    switch (rIpNum)
    {
    case 0:
        return c;
    case 1:
        return c;
    case 2:
        return c;
    case 3:
        return d;
    case 4:
        return d;
    case 5:
        return d;
    default:
        throw MechanicsException(__PRETTY_FUNCTION__, "Ip number out of range.");
    }
}


std::string NuTo::IntegrationType2D3NGauss6Ip::GetStrIdentifier() const
{
    return GetStrIdentifierStatic();
}


std::string NuTo::IntegrationType2D3NGauss6Ip::GetStrIdentifierStatic()
{
    return std::string("2D3NGAUSS6IP");
}

#ifdef ENABLE_VISUALIZE
void NuTo::IntegrationType2D3NGauss6Ip::GetVisualizationCells(unsigned int& NumVisualizationPoints, std::vector<double>& VisualizationPointLocalCoordinates, unsigned int& NumVisualizationCells, std::vector<NuTo::eCellTypes>& VisualizationCellType,
        std::vector<unsigned int>& VisualizationCellsIncidence, std::vector<unsigned int>& VisualizationCellsIP) const
{

    // only 3 integration points (1,2,3) are visualised. TODO: Voronoi decomposition + triangulation for proper visualisation

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
template void NuTo::IntegrationType2D3NGauss6Ip::serialize(boost::archive::binary_oarchive & ar, const unsigned int version);
template void NuTo::IntegrationType2D3NGauss6Ip::serialize(boost::archive::xml_oarchive & ar, const unsigned int version);
template void NuTo::IntegrationType2D3NGauss6Ip::serialize(boost::archive::text_oarchive & ar, const unsigned int version);
template void NuTo::IntegrationType2D3NGauss6Ip::serialize(boost::archive::binary_iarchive & ar, const unsigned int version);
template void NuTo::IntegrationType2D3NGauss6Ip::serialize(boost::archive::xml_iarchive & ar, const unsigned int version);
template void NuTo::IntegrationType2D3NGauss6Ip::serialize(boost::archive::text_iarchive & ar, const unsigned int version);
template<class Archive>
void NuTo::IntegrationType2D3NGauss6Ip::serialize(Archive & ar, const unsigned int version)
{
#ifdef DEBUG_SERIALIZATION
    std::cout << "start serialize IntegrationType2D3NGauss6Ip" << std::endl;
#endif
    ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(IntegrationType2D);
#ifdef DEBUG_SERIALIZATION
    std::cout << "finish serialize IntegrationType2D3NGauss6Ip" << std::endl;
#endif
}
BOOST_CLASS_EXPORT_IMPLEMENT(NuTo::IntegrationType2D3NGauss6Ip)
#endif // ENABLE_SERIALIZATION
