// $Id$

#ifdef ENABLE_VISUALIZE
#include "visualize/VisualizeEnum.h"
#endif // ENABLE_VISUALIZE

#include "mechanics/integrationtypes/IntegrationType2D4NConstVariableIp.h"

#include <assert.h>

#include <cmath>
#include <sstream>
#include <iostream>

#include "mechanics/MechanicsException.h"

//! @brief constructor
NuTo::IntegrationType2D4NConstVariableIp::IntegrationType2D4NConstVariableIp(int rNumIp)
{
    if (rNumIp<1)
        throw MechanicsException("[NuTo::IntegrationType1D2NConstVariableIp] Number of integration points must be positive.");

    mNumIp = ((int)sqrt((double)rNumIp))*((int)sqrt((double)rNumIp));
}

//! @brief returns the local coordinates of an integration point
//! @param rIpNum integration point (counting from zero)
//! @param rCoordinates (result)
void NuTo::IntegrationType2D4NConstVariableIp::GetLocalIntegrationPointCoordinates2D(int rIpNum, double rCoordinates[2])const
{
    assert(rIpNum>=0 && rIpNum<mNumIp);
	// for(int i=0; i<numIp; ++i)
    // {
    // 	std::cout << "coordinates: [" 	<< (0.5+rIpNum%(int)sqrt((double)mNumIp))/(sqrt((double)mNumIp))*2-1 <<
    // 								";" << (0.5+rIpNum/(int)sqrt((double)mNumIp))/(sqrt((double)mNumIp))*2-1 << "]" << std::endl;
    // }

    rCoordinates[0] = (0.5+rIpNum%(int)sqrt((double)mNumIp))/(sqrt((double)mNumIp))*2-1;
    rCoordinates[1] = (0.5+rIpNum/(int)sqrt((double)mNumIp))/(sqrt((double)mNumIp))*2-1;

}

//! @brief returns the total number of integration points for this integration type
//! @return number of integration points
int NuTo::IntegrationType2D4NConstVariableIp::GetNumIntegrationPoints()const
{
    return mNumIp;
}

//! @brief returns the weight of an integration point
//! @param rIpNum integration point (counting from zero)
//! @return weight of integration points
double NuTo::IntegrationType2D4NConstVariableIp::GetIntegrationPointWeight(int rIpNum)const
{
    return 4./mNumIp;
}

//! @brief returns a string with the identifier of the integration type
//! @return identifier
std::string NuTo::IntegrationType2D4NConstVariableIp::GetStrIdentifier()const
{
    std::ostringstream o;
    o << mNumIp;
    return std::string("2D4NConst"+ o.str() + "Ip");
}

#ifdef ENABLE_VISUALIZE
void NuTo::IntegrationType2D4NConstVariableIp::GetVisualizationCells(
    unsigned int& NumVisualizationPoints,
    std::vector<double>& VisualizationPointLocalCoordinates,
    unsigned int& NumVisualizationCells,
    std::vector<NuTo::eCellTypes>& VisualizationCellType,
    std::vector<unsigned int>& VisualizationCellsIncidence,
    std::vector<unsigned int>& VisualizationCellsIP) const
{
	unsigned int mNumIp1D=(unsigned int)sqrt((double)mNumIp);
	NumVisualizationPoints = (mNumIp1D+1)*(mNumIp1D+1);
    for(unsigned int i=0; i<NumVisualizationPoints; ++i)
    {
		//std::cout << "VisulizeCoordinates: [" 	<< (i%((int)sqrt((double)mNumIp)+1))/(sqrt((double)mNumIp))*2-1 <<
		//									";" << (i/((int)sqrt((double)mNumIp)+1))/(sqrt((double)mNumIp))*2-1 << "]" << std::endl;
        VisualizationPointLocalCoordinates.push_back( (i%(mNumIp1D+1))/((double)mNumIp1D)*2-1 );
        VisualizationPointLocalCoordinates.push_back( (i/(mNumIp1D+1))/((double)mNumIp1D)*2-1 );
    }
    NumVisualizationCells = mNumIp;
    unsigned int thisCell=0;
    for(unsigned int i=0; i<mNumIp1D; ++i)
    {
        for(unsigned int j=0; j<mNumIp1D; ++j)
        {
			VisualizationCellType.push_back(NuTo::eCellTypes::QUAD);
			VisualizationCellsIncidence.push_back(i*(mNumIp1D+1)+j               );
			VisualizationCellsIncidence.push_back(i*(mNumIp1D+1)+j+1             );
			VisualizationCellsIncidence.push_back(i*(mNumIp1D+1)+j+1+(mNumIp1D+1));
			VisualizationCellsIncidence.push_back(i*(mNumIp1D+1)+j  +(mNumIp1D+1));
			VisualizationCellsIP.push_back(thisCell++);
		}
    }
}
#endif // ENABLE_VISUALIZE
