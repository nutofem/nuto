#ifdef ENABLE_VISUALIZE
#include "nuto/visualize/VisualizeEnum.h"
#endif // ENABLE_VISUALIZE

#include "nuto/mechanics/integrationtypes/IntegrationType2DMod.h"

#include <assert.h>
#include <string>
#include <vector>

#include <boost/foreach.hpp>
#include <boost/utility.hpp>

#include "nuto/mechanics/MechanicsException.h"
#include "nuto/mechanics/integrationtypes/IntegrationPointBase.h"

NuTo::IntegrationType2DMod::IntegrationType2DMod()
{
    mName='\0';
}


NuTo::IntegrationType2DMod::IntegrationType2DMod(std::string rName): mName(rName)
{
}


NuTo::IntegrationType2DMod::~IntegrationType2DMod()
{
}


void NuTo::IntegrationType2DMod::GetLocalIntegrationPointCoordinates2D(int rIpNum, double rCoordinates[2]) const
{
	assert(rIpNum >= 0 && rIpNum <(int)mIpMap.size());
    boost::ptr_map<int, IntegrationPointBase>::const_iterator thisIP = boost::next(mIpMap.begin(),rIpNum);
    if (thisIP == mIpMap.end())
    	throw MechanicsException(__PRETTY_FUNCTION__, "IP " + std::to_string(rIpNum) + " does exist.");

    const std::vector<double> thisCoords(thisIP->second->GetLocalCoords());
    rCoordinates[0] = thisCoords[0];
    rCoordinates[1] = thisCoords[1];

}


unsigned int NuTo::IntegrationType2DMod::GetNumIntegrationPoints() const
{
   return mIpMap.size();
}


double NuTo::IntegrationType2DMod::GetIntegrationPointWeight(int rIpNum) const
{
	assert(rIpNum >= 0 && rIpNum < (int)mIpMap.size());
    boost::ptr_map<int, IntegrationPointBase>::const_iterator thisIP = boost::next(mIpMap.begin(), rIpNum);
    if (thisIP == mIpMap.end())
    	throw MechanicsException(__PRETTY_FUNCTION__, "IP " + std::to_string(rIpNum) + " does not exist.");

    return thisIP->second->GetWeight();
}


std::string NuTo::IntegrationType2DMod::GetStrIdentifier() const
{
    return mName;
}

#ifdef ENABLE_VISUALIZE
void NuTo::IntegrationType2DMod::GetVisualizationCells(
    unsigned int& rNumVisualizationPoints,
    std::vector<double>& rVisualizationPointLocalCoordinates,
    unsigned int& rNumVisualizationCells,
    std::vector<NuTo::eCellTypes>& rVisualizationCellType,
    std::vector<unsigned int>& rVisualizationCellsIncidence,
    std::vector<unsigned int>& rVisualizationCellsIP) const
{
	//! first initialize
	rNumVisualizationPoints=0;
	rVisualizationPointLocalCoordinates=std::vector<double>(0);
	rNumVisualizationCells=0;
	rVisualizationCellType=std::vector<NuTo::eCellTypes>(0);
	rVisualizationCellsIncidence=std::vector<unsigned int>(0);
	rVisualizationCellsIP=std::vector<unsigned int>(0);

	unsigned int numVisualizationPoints=0;
	NuTo::eCellTypes visualizationCellType=NuTo::eCellTypes::QUAD;
	std::vector<double> visualizationPointLocalCoordinates(0);
	std::vector<unsigned int> visualizationCellsIncidence(0);

	unsigned int inc=0; //< the base of the incidence for this cell
	size_t newIp=0;
	typedef boost::ptr_map< int, IntegrationPointBase > IpMap_t;
	BOOST_FOREACH( IpMap_t::const_iterator::reference it,  mIpMap)
    {
		it.second->GetVisualizationCell(numVisualizationPoints, visualizationCellType,
                visualizationPointLocalCoordinates, visualizationCellsIncidence);

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


void NuTo::IntegrationType2DMod::AddIntegrationPoint(const IntegrationPointBase &rIp)
{
	//! with this informations we build up the new IP
	IntegrationPointBase* ptrIP(0);
	ptrIP = new IntegrationPointBase(rIp);
	//! and finally insert it into the map
	int id=mIpMap.size();
	this->mIpMap.insert(id,ptrIP);
}


void NuTo::IntegrationType2DMod::DeleteIntegrationPoint(const int rIpNum)
{
	// find IP
	typedef boost::ptr_map< int, IntegrationPointBase > IpMap_t;
	IpMap_t::iterator itIP = mIpMap.find(rIpNum);
    if (itIP == this->mIpMap.end())
    {
        throw MechanicsException(__PRETTY_FUNCTION__, "Integrationpoint does not exist.");
    }
    else
    {
		// delete IP from map
		this->mIpMap.erase(rIpNum);
    }
}

