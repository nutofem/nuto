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

#include "mechanics/integrationtypes/IntegrationType2D4NLobatto16Ip.h"
#include "mechanics/integrationtypes/IntegrationType1D2NLobatto4Ip.h"


//! @brief constructor
NuTo::IntegrationType2D4NLobatto16Ip::IntegrationType2D4NLobatto16Ip()
{
    NuTo::IntegrationType1D2NLobatto4Ip Lobatto1D2N4Ip;
    Eigen::VectorXd coordinates1D2N4Ip(4);
    Eigen::VectorXd weights1D2N4Ip(4);

    // get the 1D integration point coordinates and weights
    for (int i = 0; i < 4; i++)
    {
        coordinates1D2N4Ip[i] = Lobatto1D2N4Ip.GetLocalIntegrationPointCoordinates(i)[0];
        weights1D2N4Ip[i] = Lobatto1D2N4Ip.GetIntegrationPointWeight(i);
    }

    // calculate the 2D integratration point coordinates and weights
    mWeights.resize(GetNumIntegrationPoints());
    mPts.resize(GetNumIntegrationPoints());
    int ipNum = 0;
    for (int i = 0; i < 4; i++)
        for (int j= 0; j < 4; j++)
        {
            mWeights[ipNum] = weights1D2N4Ip[i]*weights1D2N4Ip[j];
            mPts[ipNum][0] = coordinates1D2N4Ip[j];
            mPts[ipNum][1] = coordinates1D2N4Ip[i];
            ipNum++;
        }
}

//! @brief returns the local coordinates of an integration point
//! @param rIpNum integration point (counting from zero)
//! @param rCoordinates (result)
Eigen::VectorXd NuTo::IntegrationType2D4NLobatto16Ip::GetLocalIntegrationPointCoordinates(int rIpNum) const
{
    if (rIpNum>=0 && rIpNum<16)
        return mPts[rIpNum];
    else
        throw MechanicsException("[NuTo::IntegrationType2D4NLobatto16Ip::GetLocalIntegrationPointCoordinates] Ip number out of range.");
}


//! @brief returns the total number of integration points for this integration type
//! @return number of integration points
int NuTo::IntegrationType2D4NLobatto16Ip::GetNumIntegrationPoints()const
{
    return 16;
}

//! @brief returns the weight of an integration point
//! @param rIpNum integration point (counting from zero)
//! @return weight of integration points
double NuTo::IntegrationType2D4NLobatto16Ip::GetIntegrationPointWeight(int rIpNum)const
{
    if (rIpNum>=0 && rIpNum<16) return mWeights[rIpNum];
    throw MechanicsException("[NuTo::IntegrationType2D4NLobatto16Ip::GetLocalIntegrationPointCoordinates] Ip number out of range.");
}

#ifdef ENABLE_VISUALIZE
void NuTo::IntegrationType2D4NLobatto16Ip::GetVisualizationCells(
    unsigned int& NumVisualizationPoints,
    std::vector<double>& VisualizationPointLocalCoordinates,
    unsigned int& NumVisualizationCells,
    std::vector<NuTo::eCellTypes>& VisualizationCellType,
    std::vector<unsigned int>& VisualizationCellsIncidence,
    std::vector<unsigned int>& VisualizationCellsIP) const
{
    NumVisualizationPoints = 25;

    // Point 0
    VisualizationPointLocalCoordinates.push_back(-1);
    VisualizationPointLocalCoordinates.push_back(-1);

    // Point 1 (accuracy of coordinates is only important for visualization
    VisualizationPointLocalCoordinates.push_back(-0.7236);
    VisualizationPointLocalCoordinates.push_back(-1);

    // Point 2
    VisualizationPointLocalCoordinates.push_back(0.0);
    VisualizationPointLocalCoordinates.push_back(-1);

    // Point 3
    VisualizationPointLocalCoordinates.push_back(0.7236);
    VisualizationPointLocalCoordinates.push_back(-1.);

    // Point 4
    VisualizationPointLocalCoordinates.push_back(1);
    VisualizationPointLocalCoordinates.push_back(-1.);


    // Point 5
    VisualizationPointLocalCoordinates.push_back(-1);
    VisualizationPointLocalCoordinates.push_back(-0.7236);

    // Point 6
    VisualizationPointLocalCoordinates.push_back(-0.7236);
    VisualizationPointLocalCoordinates.push_back(-0.7236);

    // Point 7
    VisualizationPointLocalCoordinates.push_back(0.0);
    VisualizationPointLocalCoordinates.push_back(-0.7236);

    // Point 8
    VisualizationPointLocalCoordinates.push_back(0.7236);
    VisualizationPointLocalCoordinates.push_back(-0.7236);

    // Point 9
    VisualizationPointLocalCoordinates.push_back(1);
    VisualizationPointLocalCoordinates.push_back(-0.7236);


    // Point 10
    VisualizationPointLocalCoordinates.push_back(-1);
    VisualizationPointLocalCoordinates.push_back(-0.);

    // Point 11
    VisualizationPointLocalCoordinates.push_back(-0.7236);
    VisualizationPointLocalCoordinates.push_back(-0.);

    // Point 12
    VisualizationPointLocalCoordinates.push_back(0.0);
    VisualizationPointLocalCoordinates.push_back(-0.);

    // Point 13
    VisualizationPointLocalCoordinates.push_back(0.7236);
    VisualizationPointLocalCoordinates.push_back(-0.);

    // Point 14
    VisualizationPointLocalCoordinates.push_back(1);
    VisualizationPointLocalCoordinates.push_back(-0.);


    // Point 15
    VisualizationPointLocalCoordinates.push_back(-1);
    VisualizationPointLocalCoordinates.push_back(+0.7236);

    // Point 16
    VisualizationPointLocalCoordinates.push_back(-0.7236);
    VisualizationPointLocalCoordinates.push_back(+0.7236);

    // Point 17
    VisualizationPointLocalCoordinates.push_back(0.0);
    VisualizationPointLocalCoordinates.push_back(+0.7236);

    // Point 18
    VisualizationPointLocalCoordinates.push_back(0.7236);
    VisualizationPointLocalCoordinates.push_back(+0.7236);

    // Point 19
    VisualizationPointLocalCoordinates.push_back(1);
    VisualizationPointLocalCoordinates.push_back(+0.7236);


    // Point 20
    VisualizationPointLocalCoordinates.push_back(-1);
    VisualizationPointLocalCoordinates.push_back(+1);

    // Point 21
    VisualizationPointLocalCoordinates.push_back(-0.7236);
    VisualizationPointLocalCoordinates.push_back(+1);

    // Point 22
    VisualizationPointLocalCoordinates.push_back(0.0);
    VisualizationPointLocalCoordinates.push_back(+1);

    // Point 23
    VisualizationPointLocalCoordinates.push_back(0.7236);
    VisualizationPointLocalCoordinates.push_back(+1.);

    // Point 24
    VisualizationPointLocalCoordinates.push_back(1);
    VisualizationPointLocalCoordinates.push_back(+1.);


NumVisualizationCells = 16;

    // cell 0
    VisualizationCellType.push_back(NuTo::eCellTypes::QUAD);
    VisualizationCellsIncidence.push_back(0);
    VisualizationCellsIncidence.push_back(1);
    VisualizationCellsIncidence.push_back(6);
    VisualizationCellsIncidence.push_back(5);
    VisualizationCellsIP.push_back(0);

    // cell 1
    VisualizationCellType.push_back(NuTo::eCellTypes::QUAD);
    VisualizationCellsIncidence.push_back(1);
    VisualizationCellsIncidence.push_back(2);
    VisualizationCellsIncidence.push_back(7);
    VisualizationCellsIncidence.push_back(6);
    VisualizationCellsIP.push_back(1);

    // cell 2
    VisualizationCellType.push_back(NuTo::eCellTypes::QUAD);
    VisualizationCellsIncidence.push_back(2);
    VisualizationCellsIncidence.push_back(3);
    VisualizationCellsIncidence.push_back(8);
    VisualizationCellsIncidence.push_back(7);
    VisualizationCellsIP.push_back(2);

    // cell 3
    VisualizationCellType.push_back(NuTo::eCellTypes::QUAD);
    VisualizationCellsIncidence.push_back(3);
    VisualizationCellsIncidence.push_back(4);
    VisualizationCellsIncidence.push_back(9);
    VisualizationCellsIncidence.push_back(8);
    VisualizationCellsIP.push_back(3);

    // cell 4
    VisualizationCellType.push_back(NuTo::eCellTypes::QUAD);
    VisualizationCellsIncidence.push_back(5);
    VisualizationCellsIncidence.push_back(6);
    VisualizationCellsIncidence.push_back(11);
    VisualizationCellsIncidence.push_back(10);
    VisualizationCellsIP.push_back(4);

    // cell 5
    VisualizationCellType.push_back(NuTo::eCellTypes::QUAD);
    VisualizationCellsIncidence.push_back(6);
    VisualizationCellsIncidence.push_back(7);
    VisualizationCellsIncidence.push_back(12);
    VisualizationCellsIncidence.push_back(11);
    VisualizationCellsIP.push_back(5);

    // cell 6
    VisualizationCellType.push_back(NuTo::eCellTypes::QUAD);
    VisualizationCellsIncidence.push_back(7);
    VisualizationCellsIncidence.push_back(8);
    VisualizationCellsIncidence.push_back(13);
    VisualizationCellsIncidence.push_back(12);
    VisualizationCellsIP.push_back(6);

    // cell 7
    VisualizationCellType.push_back(NuTo::eCellTypes::QUAD);
    VisualizationCellsIncidence.push_back(8);
    VisualizationCellsIncidence.push_back(9);
    VisualizationCellsIncidence.push_back(14);
    VisualizationCellsIncidence.push_back(13);
    VisualizationCellsIP.push_back(7);

    // cell 8
    VisualizationCellType.push_back(NuTo::eCellTypes::QUAD);
    VisualizationCellsIncidence.push_back(10);
    VisualizationCellsIncidence.push_back(11);
    VisualizationCellsIncidence.push_back(16);
    VisualizationCellsIncidence.push_back(15);
    VisualizationCellsIP.push_back(8);

    // cell 9
    VisualizationCellType.push_back(NuTo::eCellTypes::QUAD);
    VisualizationCellsIncidence.push_back(11);
    VisualizationCellsIncidence.push_back(12);
    VisualizationCellsIncidence.push_back(17);
    VisualizationCellsIncidence.push_back(16);
    VisualizationCellsIP.push_back(9);

    // cell 10
    VisualizationCellType.push_back(NuTo::eCellTypes::QUAD);
    VisualizationCellsIncidence.push_back(12);
    VisualizationCellsIncidence.push_back(13);
    VisualizationCellsIncidence.push_back(18);
    VisualizationCellsIncidence.push_back(17);
    VisualizationCellsIP.push_back(10);

    // cell 11
    VisualizationCellType.push_back(NuTo::eCellTypes::QUAD);
    VisualizationCellsIncidence.push_back(13);
    VisualizationCellsIncidence.push_back(14);
    VisualizationCellsIncidence.push_back(19);
    VisualizationCellsIncidence.push_back(18);
    VisualizationCellsIP.push_back(11);

    // cell 12
    VisualizationCellType.push_back(NuTo::eCellTypes::QUAD);
    VisualizationCellsIncidence.push_back(15);
    VisualizationCellsIncidence.push_back(16);
    VisualizationCellsIncidence.push_back(21);
    VisualizationCellsIncidence.push_back(20);
    VisualizationCellsIP.push_back(12);

    // cell 13
    VisualizationCellType.push_back(NuTo::eCellTypes::QUAD);
    VisualizationCellsIncidence.push_back(16);
    VisualizationCellsIncidence.push_back(17);
    VisualizationCellsIncidence.push_back(22);
    VisualizationCellsIncidence.push_back(21);
    VisualizationCellsIP.push_back(13);

    // cell 14
    VisualizationCellType.push_back(NuTo::eCellTypes::QUAD);
    VisualizationCellsIncidence.push_back(17);
    VisualizationCellsIncidence.push_back(18);
    VisualizationCellsIncidence.push_back(23);
    VisualizationCellsIncidence.push_back(22);
    VisualizationCellsIP.push_back(14);

    // cell 15
    VisualizationCellType.push_back(NuTo::eCellTypes::QUAD);
    VisualizationCellsIncidence.push_back(18);
    VisualizationCellsIncidence.push_back(19);
    VisualizationCellsIncidence.push_back(24);
    VisualizationCellsIncidence.push_back(23);
    VisualizationCellsIP.push_back(15);

}
#endif // ENABLE_VISUALIZE

#ifdef ENABLE_SERIALIZATION
// serializes the class
template void NuTo::IntegrationType2D4NLobatto16Ip::serialize(boost::archive::binary_oarchive & ar, const unsigned int version);
template void NuTo::IntegrationType2D4NLobatto16Ip::serialize(boost::archive::xml_oarchive & ar, const unsigned int version);
template void NuTo::IntegrationType2D4NLobatto16Ip::serialize(boost::archive::text_oarchive & ar, const unsigned int version);
template void NuTo::IntegrationType2D4NLobatto16Ip::serialize(boost::archive::binary_iarchive & ar, const unsigned int version);
template void NuTo::IntegrationType2D4NLobatto16Ip::serialize(boost::archive::xml_iarchive & ar, const unsigned int version);
template void NuTo::IntegrationType2D4NLobatto16Ip::serialize(boost::archive::text_iarchive & ar, const unsigned int version);
template<class Archive>
void NuTo::IntegrationType2D4NLobatto16Ip::serialize(Archive & ar, const unsigned int version)
{
#ifdef DEBUG_SERIALIZATION
    std::cout << "start serialize IntegrationType2D4NLobatto16Ip" << std::endl;
#endif
    ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(IntegrationType2D);
#ifdef DEBUG_SERIALIZATION
    std::cout << "finish serialize IntegrationType2D4NLobatto16Ip" << std::endl;
#endif
}
BOOST_CLASS_EXPORT_IMPLEMENT(NuTo::IntegrationType2D4NLobatto16Ip)
#endif // ENABLE_SERIALIZATION
