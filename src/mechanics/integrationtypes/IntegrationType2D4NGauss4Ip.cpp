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

#include "mechanics/integrationtypes/IntegrationType2D4NGauss4Ip.h"
#include <assert.h>


//! @brief constructor
NuTo::IntegrationType2D4NGauss4Ip::IntegrationType2D4NGauss4Ip()
{
}

//! @brief returns the local coordinates of an integration point
//! @param rIpNum integration point (counting from zero)
//! @param rCoordinates (result)
Eigen::VectorXd NuTo::IntegrationType2D4NGauss4Ip::GetLocalIntegrationPointCoordinates(int rIpNum) const
{
    assert(rIpNum>=0 && rIpNum<4);
    switch (rIpNum)
    {
    case 0 : return Eigen::Vector2d({-0.577350269189626, -0.577350269189626});
    case 1 : return Eigen::Vector2d({+0.577350269189626, -0.577350269189626});
    case 2 : return Eigen::Vector2d({+0.577350269189626, +0.577350269189626});
    case 3 : return Eigen::Vector2d({-0.577350269189626, +0.577350269189626});
    default:
        throw MechanicsException("[NuTo::IntegrationType2D4NGauss4Ip::GetLocalIntegrationPointCoordinates] Ip number out of range.");
    }
}


//! @brief returns the total number of integration points for this integration type
//! @return number of integration points
int NuTo::IntegrationType2D4NGauss4Ip::GetNumIntegrationPoints()const
{
    return 4;
}

//! @brief returns the weight of an integration point
//! @param rIpNum integration point (counting from zero)
//! @return weight of integration points
double NuTo::IntegrationType2D4NGauss4Ip::GetIntegrationPointWeight(int rIpNum)const
{
    return 1;
}

#ifdef ENABLE_VISUALIZE
void NuTo::IntegrationType2D4NGauss4Ip::GetVisualizationCells(
    unsigned int& NumVisualizationPoints,
    std::vector<double>& VisualizationPointLocalCoordinates,
    unsigned int& NumVisualizationCells,
    std::vector<NuTo::eCellTypes>& VisualizationCellType,
    std::vector<unsigned int>& VisualizationCellsIncidence,
    std::vector<unsigned int>& VisualizationCellsIP) const
{
    NumVisualizationPoints = 9;

    // Point 0
    VisualizationPointLocalCoordinates.push_back(-1);
    VisualizationPointLocalCoordinates.push_back(-1);

    // Point 1
    VisualizationPointLocalCoordinates.push_back(0);
    VisualizationPointLocalCoordinates.push_back(-1);

    // Point 2
    VisualizationPointLocalCoordinates.push_back(1);
    VisualizationPointLocalCoordinates.push_back(-1);

    // Point 3
    VisualizationPointLocalCoordinates.push_back(-1);
    VisualizationPointLocalCoordinates.push_back(0);

    // Point 4
    VisualizationPointLocalCoordinates.push_back(0);
    VisualizationPointLocalCoordinates.push_back(0);

    // Point 5
    VisualizationPointLocalCoordinates.push_back(1);
    VisualizationPointLocalCoordinates.push_back(0);

    // Point 6
    VisualizationPointLocalCoordinates.push_back(-1);
    VisualizationPointLocalCoordinates.push_back(+1);

    // Point 7
    VisualizationPointLocalCoordinates.push_back(0);
    VisualizationPointLocalCoordinates.push_back(+1);

    // Point 8
    VisualizationPointLocalCoordinates.push_back(1);
    VisualizationPointLocalCoordinates.push_back(+1);

    NumVisualizationCells = 4;

    // cell 0
    VisualizationCellType.push_back(NuTo::eCellTypes::QUAD);
    VisualizationCellsIncidence.push_back(0);
    VisualizationCellsIncidence.push_back(1);
    VisualizationCellsIncidence.push_back(4);
    VisualizationCellsIncidence.push_back(3);
    VisualizationCellsIP.push_back(0);

    // cell 1
    VisualizationCellType.push_back(NuTo::eCellTypes::QUAD);
    VisualizationCellsIncidence.push_back(1);
    VisualizationCellsIncidence.push_back(2);
    VisualizationCellsIncidence.push_back(5);
    VisualizationCellsIncidence.push_back(4);
    VisualizationCellsIP.push_back(1);

    // cell 2
    VisualizationCellType.push_back(NuTo::eCellTypes::QUAD);
    VisualizationCellsIncidence.push_back(4);
    VisualizationCellsIncidence.push_back(5);
    VisualizationCellsIncidence.push_back(8);
    VisualizationCellsIncidence.push_back(7);
    VisualizationCellsIP.push_back(2);

    // cell 3
    VisualizationCellType.push_back(NuTo::eCellTypes::QUAD);
    VisualizationCellsIncidence.push_back(3);
    VisualizationCellsIncidence.push_back(4);
    VisualizationCellsIncidence.push_back(7);
    VisualizationCellsIncidence.push_back(6);
    VisualizationCellsIP.push_back(3);
}
#endif // ENABLE_VISUALIZE

#ifdef ENABLE_SERIALIZATION
// serializes the class
template void NuTo::IntegrationType2D4NGauss4Ip::serialize(boost::archive::binary_oarchive & ar, const unsigned int version);
template void NuTo::IntegrationType2D4NGauss4Ip::serialize(boost::archive::xml_oarchive & ar, const unsigned int version);
template void NuTo::IntegrationType2D4NGauss4Ip::serialize(boost::archive::text_oarchive & ar, const unsigned int version);
template void NuTo::IntegrationType2D4NGauss4Ip::serialize(boost::archive::binary_iarchive & ar, const unsigned int version);
template void NuTo::IntegrationType2D4NGauss4Ip::serialize(boost::archive::xml_iarchive & ar, const unsigned int version);
template void NuTo::IntegrationType2D4NGauss4Ip::serialize(boost::archive::text_iarchive & ar, const unsigned int version);
template<class Archive>
void NuTo::IntegrationType2D4NGauss4Ip::serialize(Archive & ar, const unsigned int version)
{
#ifdef DEBUG_SERIALIZATION
    std::cout << "start serialize IntegrationType2D4NGauss4Ip" << std::endl;
#endif
    ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(IntegrationType2D);
#ifdef DEBUG_SERIALIZATION
    std::cout << "finish serialize IntegrationType2D4NGauss4Ip" << std::endl;
#endif
}
BOOST_CLASS_EXPORT_IMPLEMENT(NuTo::IntegrationType2D4NGauss4Ip)
#endif // ENABLE_SERIALIZATION
