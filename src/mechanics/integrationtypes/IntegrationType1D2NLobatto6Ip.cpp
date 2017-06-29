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

#include "mechanics/integrationtypes/IntegrationType1D2NLobatto6Ip.h"

// constructor
NuTo::IntegrationType1D2NLobatto6Ip::IntegrationType1D2NLobatto6Ip():
    iPts{-1., -0.76505532392946469, -0.28523151648064509, 0.28523151648064509, +0.76505532392946469, +1.},
    weights{0.066666666666666666, 0.554858377035486353, 0.378474956297846980, 0.378474956297846980, 0.554858377035486353, 0.066666666666666666}
{
}

//! @brief returns the local coordinates of an integration point
//! @param rIpNum integration point (counting from zero)
//! @param rCoordinates (result)
Eigen::VectorXd NuTo::IntegrationType1D2NLobatto6Ip::GetLocalIntegrationPointCoordinates(int rIpNum) const
{
    if(rIpNum >= 0 && rIpNum < 6)
        return Eigen::Matrix<double, 1, 1>::Constant(iPts[rIpNum]);
    else
        throw Exception("[NuTo::IntegrationType1D2NLobatto6Ip::GetLocalIntegrationPointCoordinates] Ip number out of range.");
}


//! @brief returns the total number of integration points for this integration type
//! @return number of integration points
int NuTo::IntegrationType1D2NLobatto6Ip::GetNumIntegrationPoints()const
{
    return 6;
}

//! @brief returns the weight of an integration point
//! @param rIpNum integration point (counting from zero)
//! @return weight of integration points
double NuTo::IntegrationType1D2NLobatto6Ip::GetIntegrationPointWeight(int rIpNum)const
{
    if(rIpNum >= 0 && rIpNum < 6) return weights[rIpNum];
    throw Exception("[NuTo::IntegrationType1D2NLobatto6Ip::GetIntegrationPointWeight] Ip number out of range.");
}

#ifdef ENABLE_VISUALIZE
void NuTo::IntegrationType1D2NLobatto6Ip::GetVisualizationCells(
    unsigned int& NumVisualizationPoints,
    std::vector<double>& VisualizationPointLocalCoordinates,
    unsigned int& NumVisualizationCells,
    std::vector<NuTo::eCellTypes>& VisualizationCellType,
    std::vector<unsigned int>& VisualizationCellsIncidence,
    std::vector<unsigned int>& VisualizationCellsIP) const
{
    NumVisualizationPoints = 7;
    VisualizationPointLocalCoordinates.push_back(-1.);
    VisualizationPointLocalCoordinates.push_back(-0.88252);
    VisualizationPointLocalCoordinates.push_back(-0.52514);
    VisualizationPointLocalCoordinates.push_back( 0);
    VisualizationPointLocalCoordinates.push_back( 0.88252);
    VisualizationPointLocalCoordinates.push_back(0.52514);
    VisualizationPointLocalCoordinates.push_back( 1.);
    NumVisualizationCells = 6;
    VisualizationCellType.push_back(NuTo::eCellTypes::LINE);
    VisualizationCellType.push_back(NuTo::eCellTypes::LINE);
    VisualizationCellType.push_back(NuTo::eCellTypes::LINE);
    VisualizationCellType.push_back(NuTo::eCellTypes::LINE);
    VisualizationCellType.push_back(NuTo::eCellTypes::LINE);
    VisualizationCellType.push_back(NuTo::eCellTypes::LINE);

    VisualizationCellsIncidence.push_back(0);
    VisualizationCellsIncidence.push_back(1);

    VisualizationCellsIncidence.push_back(1);
    VisualizationCellsIncidence.push_back(2);

    VisualizationCellsIncidence.push_back(2);
    VisualizationCellsIncidence.push_back(3);

    VisualizationCellsIncidence.push_back(3);
    VisualizationCellsIncidence.push_back(4);

    VisualizationCellsIncidence.push_back(4);
    VisualizationCellsIncidence.push_back(5);

    VisualizationCellsIncidence.push_back(5);
    VisualizationCellsIncidence.push_back(6);

    VisualizationCellsIP.push_back(0);
    VisualizationCellsIP.push_back(1);
    VisualizationCellsIP.push_back(2);
    VisualizationCellsIP.push_back(3);
    VisualizationCellsIP.push_back(4);
    VisualizationCellsIP.push_back(5);
}
#endif // ENABLE_VISUALIZE

#ifdef ENABLE_SERIALIZATION
BOOST_CLASS_EXPORT_IMPLEMENT(NuTo::IntegrationType1D2NLobatto6Ip)
#endif
