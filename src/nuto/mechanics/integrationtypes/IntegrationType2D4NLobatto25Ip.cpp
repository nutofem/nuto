// $Id: IntegrationType2D4NLobatto25Ip.cpp 331 2010-10-06 09:32:11Z arnold2 $

#ifdef ENABLE_SERIALIZATION
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/xml_oarchive.hpp>
#include <boost/archive/xml_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#endif //ENABLE_SERIALIZATION

#include "nuto/mechanics/integrationtypes/IntegrationType2D4NLobatto25Ip.h"
#include <assert.h>

//! @brief constructor
NuTo::IntegrationType2D4NLobatto25Ip::IntegrationType2D4NLobatto25Ip()
{
}


//! @brief returns the local coordinates of an integration point
//! @param rIpNum integration point (counting from zero)
//! @param rCoordinates (result)
void NuTo::IntegrationType2D4NLobatto25Ip::GetLocalIntegrationPointCoordinates2D(int rIpNum, double rCoordinates[2])const
{
    assert(rIpNum>=0 && rIpNum<25);
    switch (rIpNum)
    {
    case 0 :
        rCoordinates[0] = -1.0;
        rCoordinates[1] = -1.0;
        break;
    case 1 :
        rCoordinates[0] = -0.654653670707977087;
        rCoordinates[1] = -1.0;
        break;
    case 2 :
        rCoordinates[0] = +0.0;
        rCoordinates[1] = -1.0;
        break;
    case 3 :
        rCoordinates[0] = +0.654653670707977087;
        rCoordinates[1] = -1.0;
        break;
    case 4 :
        rCoordinates[0] = +1.0;
        rCoordinates[1] = -1;
        break;

    case 5 :
        rCoordinates[0] = -1.0;
        rCoordinates[1] = -0.654653670707977087;
        break;
    case 6 :
        rCoordinates[0] = -0.654653670707977087;
        rCoordinates[1] = -0.654653670707977087;
        break;
    case 7 :
        rCoordinates[0] = +0.0;
        rCoordinates[1] = -0.654653670707977087;
        break;
    case 8 :
        rCoordinates[0] = +0.654653670707977087;
        rCoordinates[1] = -0.654653670707977087;
        break;
    case 9 :
        rCoordinates[0] = +1.0;
        rCoordinates[1] = -0.654653670707977087;
        break;

    case 10 :
        rCoordinates[0] = -1.0;
        rCoordinates[1] = 0.0;
        break;
    case 11 :
        rCoordinates[0] = -0.654653670707977087;
        rCoordinates[1] = 0.0;
        break;
    case 12 :
        rCoordinates[0] = +0.0;
        rCoordinates[1] = 0.0;
        break;
    case 13 :
        rCoordinates[0] = +0.654653670707977087;
        rCoordinates[1] = 0.0;
        break;
    case 14 :
        rCoordinates[0] = +1.0;
        rCoordinates[1] = 0.0;
        break;

    case 15 :
        rCoordinates[0] = -1.0;
        rCoordinates[1] = +0.654653670707977087;
        break;
    case 16 :
        rCoordinates[0] = -0.654653670707977087;
        rCoordinates[1] = +0.654653670707977087;
        break;
    case 17 :
        rCoordinates[0] = +0.0;
        rCoordinates[1] = +0.654653670707977087;
        break;
    case 18 :
        rCoordinates[0] = +0.654653670707977087;
        rCoordinates[1] = +0.654653670707977087;
        break;
    case 19 :
        rCoordinates[0] = +1.0;
        rCoordinates[1] = +0.654653670707977087;
        break;

    case 20 :
        rCoordinates[0] = -1.0;
        rCoordinates[1] = 1.0;
        break;
    case 21 :
        rCoordinates[0] = -0.654653670707977087;
        rCoordinates[1] = 1.0;
        break;
    case 22 :
        rCoordinates[0] = +0.0;
        rCoordinates[1] = 1.0;
        break;
    case 23 :
        rCoordinates[0] = +0.654653670707977087;
        rCoordinates[1] = 1.0;
        break;
    case 24 :
        rCoordinates[0] = +1.0;
        rCoordinates[1] = 1.0;
        break;
    default:
        throw MechanicsException("[NuTo::IntegrationType2D4NLobatto25Ip::GetLocalIntegrationPointCoordinates] Ip number out of range.");
    }
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
    switch (rIpNum)
    {
    case 0 :
        return 0.01;
        break;
    case 1 :
    	return 0.05444444444444444;
        break;
    case 2 :
    	return 0.07111111111111111;
        break;
    case 3 :
    	return 0.05444444444444444;
        break;
    case 4 :
    	return 0.01;
        break;
    case 5 :
    	return 0.05444444444444444;
        break;
    case 6 :
    	return 0.296419753086419713;
        break;
    case 7 :
    	return 0.387160493827160501;
        break;
    case 8 :
    	return 0.296419753086419713;
        break;
    case 9 :
    	return 0.054444444444444441;
        break;
    case 10 :
    	return 0.071111111111111111;
        break;
    case 11 :
    	return 0.387160493827160501;
        break;
    case 12 :
    	return 0.505679012345679024;
        break;
    case 13 :
    	return 0.387160493827160501;
        break;
    case 14 :
    	return 0.071111111111111111;
        break;
    case 15 :
    	return 0.05444444444444444;
        break;
    case 16 :
    	return 0.296419753086419713;
        break;
    case 17 :
    	return 0.387160493827160501;
        break;
    case 18 :
    	return 0.296419753086419713;
        break;
    case 19 :
    	return 0.054444444444444441;
        break;
    case 20 :
    	return 0.01;
        break;
    case 21 :
    	return 0.05444444444444444;
        break;
    case 22 :
    	return 0.071111111111111111;
        break;
    case 23 :
    	return 0.05444444444444444;
        break;
    case 24 :
    	return 0.01;
        break;
    default:
        throw MechanicsException("[NuTo::IntegrationType2D4NLobatto25Ip::GetLocalIntegrationPointCoordinates] Ip number out of range.");
    }
}

//! @brief returns a string with the identifier of the integration type
//! @return identifier
std::string NuTo::IntegrationType2D4NLobatto25Ip::GetStrIdentifier()const
{
    return GetStrIdentifierStatic();
}

//! @brief returns a string with the identifier of the integration type
//! @return identifier
std::string NuTo::IntegrationType2D4NLobatto25Ip::GetStrIdentifierStatic()
{
    return std::string("2D4NLOBATTO25IP");
}

#ifdef ENABLE_VISUALIZE
void NuTo::IntegrationType2D4NLobatto25Ip::GetVisualizationCells(
    unsigned int& NumVisualizationPoints,
    std::vector<double>& VisualizationPointLocalCoordinates,
    unsigned int& NumVisualizationCells,
    std::vector<NuTo::CellBase::eCellTypes>& VisualizationCellType,
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
    	    VisualizationCellType.push_back(NuTo::CellBase::QUAD);
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
