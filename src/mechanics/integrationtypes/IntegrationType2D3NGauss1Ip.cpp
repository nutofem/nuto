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

#include "mechanics/integrationtypes/IntegrationType2D3NGauss1Ip.h"


//! @brief constructor
NuTo::IntegrationType2D3NGauss1Ip::IntegrationType2D3NGauss1Ip()
{
}

//! @brief returns the local coordinates of an integration point
//! @param rIpNum integration point (counting from zero)
//! @param rCoordinates (result)
Eigen::VectorXd NuTo::IntegrationType2D3NGauss1Ip::GetLocalIntegrationPointCoordinates(int rIpNum) const
{
    if (rIpNum!=0)
        throw Exception("[NuTo::IntegrationType2D3NGauss1Ip::GetLocalIntegrationPointCoordinates] Ip number out of range.");

    return Eigen::Vector2d({1./3., 1./3.});
}


//! @brief returns the total number of integration points for this integration type
//! @return number of integration points
int NuTo::IntegrationType2D3NGauss1Ip::GetNumIntegrationPoints()const
{
    return 1;
}

//! @brief returns the weight of an integration point
//! @param rIpNum integration point (counting from zero)
//! @return weight of integration points
double NuTo::IntegrationType2D3NGauss1Ip::GetIntegrationPointWeight(int rIpNum)const
{
    return 0.5;
}

#ifdef ENABLE_VISUALIZE
void NuTo::IntegrationType2D3NGauss1Ip::GetVisualizationCells(
    unsigned int& NumVisualizationPoints,
    std::vector<double>& VisualizationPointLocalCoordinates,
    unsigned int& NumVisualizationCells,
    std::vector<NuTo::eCellTypes>& VisualizationCellType,
    std::vector<unsigned int>& VisualizationCellsIncidence,
    std::vector<unsigned int>& VisualizationCellsIP) const
{
    NumVisualizationPoints = 3;

    // Point 0
    VisualizationPointLocalCoordinates.push_back(0);
    VisualizationPointLocalCoordinates.push_back(0);

    // Point 1
    VisualizationPointLocalCoordinates.push_back(1.0);
    VisualizationPointLocalCoordinates.push_back(0);

    // Point 2
    VisualizationPointLocalCoordinates.push_back(0.);
    VisualizationPointLocalCoordinates.push_back(1.0);

    NumVisualizationCells = 1;

    // cell 0
    VisualizationCellType.push_back(NuTo::eCellTypes::TRIANGLE);
    VisualizationCellsIncidence.push_back(0);
    VisualizationCellsIncidence.push_back(1);
    VisualizationCellsIncidence.push_back(2);
    VisualizationCellsIP.push_back(0);
}

#endif // ENABLE_VISUALIZE

#ifdef ENABLE_SERIALIZATION
// serializes the class
template void NuTo::IntegrationType2D3NGauss1Ip::serialize(boost::archive::binary_oarchive & ar, const unsigned int version);
template void NuTo::IntegrationType2D3NGauss1Ip::serialize(boost::archive::xml_oarchive & ar, const unsigned int version);
template void NuTo::IntegrationType2D3NGauss1Ip::serialize(boost::archive::text_oarchive & ar, const unsigned int version);
template void NuTo::IntegrationType2D3NGauss1Ip::serialize(boost::archive::binary_iarchive & ar, const unsigned int version);
template void NuTo::IntegrationType2D3NGauss1Ip::serialize(boost::archive::xml_iarchive & ar, const unsigned int version);
template void NuTo::IntegrationType2D3NGauss1Ip::serialize(boost::archive::text_iarchive & ar, const unsigned int version);
template<class Archive>
void NuTo::IntegrationType2D3NGauss1Ip::serialize(Archive & ar, const unsigned int version)
{
#ifdef DEBUG_SERIALIZATION
    std::cout << "start serialize IntegrationType2D3NGauss1Ip" << std::endl;
#endif
    ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(IntegrationType2D);
#ifdef DEBUG_SERIALIZATION
    std::cout << "finish serialize IntegrationType2D3NGauss1Ip" << std::endl;
#endif
}
BOOST_CLASS_EXPORT_IMPLEMENT(NuTo::IntegrationType2D3NGauss1Ip)
#endif // ENABLE_SERIALIZATION

