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

#include "mechanics/integrationtypes/IntegrationType2D4NLobatto25Ip.h"
#include "mechanics/integrationtypes/IntegrationType1D2NLobatto.h"

//! @brief constructor
NuTo::IntegrationType2D4NLobatto25Ip::IntegrationType2D4NLobatto25Ip()
{
    NuTo::IntegrationType1D2NLobatto Lobatto1D2N5Ip(5);
    Eigen::VectorXd coordinates1D2N5Ip(5);
	Eigen::VectorXd weights1D2N5Ip(5);

    // get the 1D integration point coordinates and weights
    for (int i = 0; i < 5; i++)
    {
		coordinates1D2N5Ip[i] = Lobatto1D2N5Ip.GetLocalIntegrationPointCoordinates(i)[0];
        weights1D2N5Ip[i] = Lobatto1D2N5Ip.GetIntegrationPointWeight(i);
    }

    // calculate the 2D integratration point coordinates and weights
	mWeights.resize(GetNumIntegrationPoints());
	mPts.resize(GetNumIntegrationPoints());
    int ipNum = 0;
    for (int i = 0; i < 5; i++)
        for (int j= 0; j < 5; j++)
        {
            mWeights[ipNum] = weights1D2N5Ip[i]*weights1D2N5Ip[j];
            mPts[ipNum][0] = coordinates1D2N5Ip[j];
            mPts[ipNum][1] = coordinates1D2N5Ip[i];
            ipNum++;
        }
}


//! @brief returns the local coordinates of an integration point
//! @param rIpNum integration point (counting from zero)
//! @param rCoordinates (result)
Eigen::VectorXd NuTo::IntegrationType2D4NLobatto25Ip::GetLocalIntegrationPointCoordinates(int rIpNum) const
{
    if (rIpNum>=0 && rIpNum<25)
		return mPts[rIpNum];
    else
        throw MechanicsException("[NuTo::IntegrationType2D4NLobatto25Ip::GetLocalIntegrationPointCoordinates] Ip number out of range.");
}


//! @brief returns the total number of integration points for this integration type
//! @return number of integration points
int NuTo::IntegrationType2D4NLobatto25Ip::GetNumIntegrationPoints()const
{
    return 25;
}

//! @brief returns the weight of an integration point
//! @param rIpNum integration point (counting from zero)
//! @return weight of integration points
double NuTo::IntegrationType2D4NLobatto25Ip::GetIntegrationPointWeight(int rIpNum)const
{
    if (rIpNum>=0 && rIpNum<25) return mWeights[rIpNum];
    throw MechanicsException("[NuTo::IntegrationType2D4NLobatto25Ip::GetIntegrationPointWeight] Ip number out of range.");
}


#ifdef ENABLE_VISUALIZE
void NuTo::IntegrationType2D4NLobatto25Ip::GetVisualizationCells(
    unsigned int& NumVisualizationPoints,
    std::vector<double>& VisualizationPointLocalCoordinates,
    unsigned int& NumVisualizationCells,
    std::vector<NuTo::eCellTypes>& VisualizationCellType,
    std::vector<unsigned int>& VisualizationCellsIncidence,
    std::vector<unsigned int>& VisualizationCellsIP) const
{
    NumVisualizationPoints = 36;

    for (int county=0; county<6; county++)
    {
    	double y;
    	switch(county)
    	{
    	case 0:
    		y=-1;
    		break;
    	case 1:
    		y=-0.83;
    		break;
    	case 2:
    		y=-0.33;
    		break;
    	case 3:
    		y = 0.33;
    		break;
    	case 4:
    		y = 0.83;
    		break;
    	case 5:
    		y = 1;
    		break;
    	}
    	for (int countx=0; countx<6; countx++)
    	{
        	double x;
        	switch(countx)
        	{
        	case 0:
        		x=-1;
        		break;
        	case 1:
        		x=-0.83;
        		break;
        	case 2:
        		x=-0.33;
        		break;
        	case 3:
        		x = 0.33;
        		break;
        	case 4:
        		x = 0.83;
        		break;
        	case 5:
        		x = 1;
        		break;
        	}
            VisualizationPointLocalCoordinates.push_back(x);
            VisualizationPointLocalCoordinates.push_back(y);
    	}
    }

    NumVisualizationCells = 25;

    int theCell(0);
    for (int county=0; county<5; county++)
    {
    	for (int countx=0; countx<5; countx++)
    	{
    	    VisualizationCellType.push_back(NuTo::eCellTypes::QUAD);
    	    VisualizationCellsIncidence.push_back(county*6+countx);
    	    VisualizationCellsIncidence.push_back(county*6+countx+1);
    	    VisualizationCellsIncidence.push_back((county+1)*6+countx+1);
    	    VisualizationCellsIncidence.push_back((county+1)*6+countx);
    	    VisualizationCellsIP.push_back(theCell);
    	    theCell++;
    	}
    }
}
#endif // ENABLE_VISUALIZE

#ifdef ENABLE_SERIALIZATION
// serializes the class
template void NuTo::IntegrationType2D4NLobatto25Ip::serialize(boost::archive::binary_oarchive & ar, const unsigned int version);
template void NuTo::IntegrationType2D4NLobatto25Ip::serialize(boost::archive::xml_oarchive & ar, const unsigned int version);
template void NuTo::IntegrationType2D4NLobatto25Ip::serialize(boost::archive::text_oarchive & ar, const unsigned int version);
template void NuTo::IntegrationType2D4NLobatto25Ip::serialize(boost::archive::binary_iarchive & ar, const unsigned int version);
template void NuTo::IntegrationType2D4NLobatto25Ip::serialize(boost::archive::xml_iarchive & ar, const unsigned int version);
template void NuTo::IntegrationType2D4NLobatto25Ip::serialize(boost::archive::text_iarchive & ar, const unsigned int version);
template<class Archive>
void NuTo::IntegrationType2D4NLobatto25Ip::serialize(Archive & ar, const unsigned int version)
{
#ifdef DEBUG_SERIALIZATION
    std::cout << "start serialize IntegrationType2D4NLobatto25Ip" << std::endl;
#endif
    ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(IntegrationType2D);
#ifdef DEBUG_SERIALIZATION
    std::cout << "finish serialize IntegrationType2D4NLobatto25Ip" << std::endl;
#endif
}
BOOST_CLASS_EXPORT_IMPLEMENT(NuTo::IntegrationType2D4NLobatto25Ip)
#endif // ENABLE_SERIALIZATION
