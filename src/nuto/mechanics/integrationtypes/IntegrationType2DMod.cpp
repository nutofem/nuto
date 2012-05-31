// $Id$

#include "nuto/mechanics/integrationtypes/IntegrationType2DMod.h"

#include <assert.h>

#include <sstream>
#include <iostream>
#include <vector>

#include <boost/foreach.hpp>
#include <boost/utility.hpp>

#include "nuto/mechanics/MechanicsException.h"
#include "nuto/mechanics/integrationtypes/IntegrationPointBase.h"

//! @brief constructor
NuTo::IntegrationType2DMod::IntegrationType2DMod()
{
    mName='\0';
}

//! @brief constructor
NuTo::IntegrationType2DMod::IntegrationType2DMod(std::string rName):
	mName(rName)
{
}

//! @brief destructor
NuTo::IntegrationType2DMod::~IntegrationType2DMod()
{
}

//! @brief returns the local coordinates of an integration point
//! @param rIpNum integration point (counting from zero)
//! @param rCoordinates (result)
void NuTo::IntegrationType2DMod::GetLocalIntegrationPointCoordinates2D(int rIpNum, double rCoordinates[2]) const
{
	assert(rIpNum>=0 && rIpNum<(int)mIpMap.size());
    boost::ptr_map< int, IntegrationPointBase >::const_iterator thisIP = boost::next(mIpMap.begin(),rIpNum);
    if (thisIP==mIpMap.end())
    {
    	std::stringstream message;
    	message << "[NuTo::IntegrationType2DMod::GetLocalIntegrationPointCoordinates2D] IP " << rIpNum << " does not exist." << std::endl;
    	throw MechanicsException(message.str());
    }

    const std::vector<double> thisCoords(thisIP->second->GetLocalCoords());
//	assert(rIpNum>=0 && rIpNum<mNumIp);
//    boost::ptr_map< int, IntegrationPointBase >::const_iterator thisIP = mIpMap.find(rIpNum);
//    if (thisIP==mIpMap.end())
//    {
//    	std::stringstream message;
//    	message << "[NuTo::IntegrationType2DMod::GetLocalIntegrationPointCoordinates2D] IP " << rIpNum << " does not exist." << std::endl;
//    	throw MechanicsException(message.str());
//    }
//
//    const std::vector<double> thisCoords(thisIP->second->GetLocalCoords());
    rCoordinates[0] = thisCoords[0];
    rCoordinates[1] = thisCoords[1];

}

//! @brief returns the total number of integration points for this integration type
//! @return number of integration points
int NuTo::IntegrationType2DMod::GetNumIntegrationPoints()const
{
   return mIpMap.size();
}

//! @brief returns the weight of an integration point
//! @param rIpNum integration point (counting from zero)
//! @return weight of integration points
double NuTo::IntegrationType2DMod::GetIntegrationPointWeight(int rIpNum)const
{
	assert(rIpNum>=0 && rIpNum<(int)mIpMap.size());
    boost::ptr_map< int, IntegrationPointBase >::const_iterator thisIP = boost::next(mIpMap.begin(),rIpNum);
    if (thisIP==mIpMap.end())
    {
    	std::stringstream message;
    	message << "[NuTo::IntegrationType2DMod::GetIntegrationPointWeight] IP " << rIpNum << " does not exist." << std::endl;
    	throw MechanicsException(message.str());
    }

    return thisIP->second->GetWeight();
}

//! @brief returns a string with the identifier of the integration type
//! @return identifier
std::string NuTo::IntegrationType2DMod::GetStrIdentifier()const
{
    return mName;
}

#ifdef ENABLE_VISUALIZE
void NuTo::IntegrationType2DMod::GetVisualizationCells(
    unsigned int& rNumVisualizationPoints,
    std::vector<double>& rVisualizationPointLocalCoordinates,
    unsigned int& rNumVisualizationCells,
    std::vector<NuTo::CellBase::eCellTypes>& rVisualizationCellType,
    std::vector<unsigned int>& rVisualizationCellsIncidence,
    std::vector<unsigned int>& rVisualizationCellsIP) const
{
	//! first initialize
	rNumVisualizationPoints=0;
	rVisualizationPointLocalCoordinates=std::vector<double>(0);
	rNumVisualizationCells=0;
	rVisualizationCellType=std::vector<NuTo::CellBase::eCellTypes>(0);
	rVisualizationCellsIncidence=std::vector<unsigned int>(0);
	rVisualizationCellsIP=std::vector<unsigned int>(0);

	unsigned int numVisualizationPoints=0;
	NuTo::CellBase::eCellTypes visualizationCellType=NuTo::CellBase::eCellTypes::QUAD;
	std::vector<double> visualizationPointLocalCoordinates(0);
	std::vector<unsigned int> visualizationCellsIncidence(0);

	unsigned int inc=0; //< the base of the incidence for this cell
	size_t newIp=0;
	typedef boost::ptr_map< int, IntegrationPointBase > IpMap_t;
	BOOST_FOREACH( IpMap_t::const_iterator::reference it,  mIpMap){
		it.second->GetVisualizationCell( numVisualizationPoints, visualizationCellType,visualizationPointLocalCoordinates, visualizationCellsIncidence );

		//! NumVisualizationPoints
		rNumVisualizationPoints+=numVisualizationPoints;

		//! VisualizationCellsIncidence
		BOOST_FOREACH( unsigned int u, visualizationCellsIncidence )
			rVisualizationCellsIncidence.push_back(inc+u);
		inc+=numVisualizationPoints;

		//! VisualizationPointLocalCoordinates
		BOOST_FOREACH( double d, visualizationPointLocalCoordinates )
			rVisualizationPointLocalCoordinates.push_back(d);

		//! VisualizationCellsIP
		rVisualizationCellsIP.push_back(newIp++);

		//! VisualizationCellType
		rVisualizationCellType.push_back(visualizationCellType);
	}

	//! NumVisualizationCells
	rNumVisualizationCells=rVisualizationCellsIP.size();
}
#endif // ENABLE_VISUALIZE

    //! @brief adds a new integration point
    //! @param rIp (Input) integration point
void NuTo::IntegrationType2DMod::AddIntegrationPoint(const IntegrationPointBase &rIp){
	//! with this informations we build up the new IP
	IntegrationPointBase* ptrIP(0);
	ptrIP = new IntegrationPointBase(rIp);
	//! and finally insert it into the map
	int id=mIpMap.size();
	this->mIpMap.insert(id,ptrIP);
}

//! @brief deletes an integration point
//! @param rIpNum (Input) integration point (counting from zero)
void NuTo::IntegrationType2DMod::DeleteIntegrationPoint(const int rIpNum)
{
	// find IP
	typedef boost::ptr_map< int, IntegrationPointBase > IpMap_t;
	IpMap_t::iterator itIP = mIpMap.find(rIpNum);
    if (itIP == this->mIpMap.end())
    {
        throw MechanicsException("[NuTo::IntegrationType2DMod::DeleteIntegrationPoint] Integrationpoint does not exist.");
    }
    else
    {
		// delete IP from map
		this->mIpMap.erase(rIpNum);
    }

}

