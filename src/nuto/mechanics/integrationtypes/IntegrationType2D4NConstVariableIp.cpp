#ifdef ENABLE_VISUALIZE
#include "nuto/visualize/VisualizeEnum.h"
#endif // ENABLE_VISUALIZE

#include "nuto/mechanics/integrationtypes/IntegrationType2D4NConstVariableIp.h"

#include <assert.h>
#include <string>
#include <cmath>

#include "nuto/mechanics/MechanicsException.h"

NuTo::IntegrationType2D4NConstVariableIp::IntegrationType2D4NConstVariableIp(int rNumIp)
{
    if (rNumIp<1)
        throw MechanicsException(__PRETTY_FUNCTION__, "Number of integration points must be positive.");

    mNumIp = ((int)sqrt((double)rNumIp))*((int)sqrt((double)rNumIp));
}


void NuTo::IntegrationType2D4NConstVariableIp::GetLocalIntegrationPointCoordinates2D(int rIpNum, double rCoordinates[2]) const
{
    assert(rIpNum>=0 && rIpNum<mNumIp);
    rCoordinates[0] = (0.5+rIpNum%(int)sqrt((double)mNumIp))/(sqrt((double)mNumIp))*2-1;
    rCoordinates[1] = (0.5+rIpNum/(int)sqrt((double)mNumIp))/(sqrt((double)mNumIp))*2-1;

}


unsigned int NuTo::IntegrationType2D4NConstVariableIp::GetNumIntegrationPoints() const
{
    return mNumIp;
}


double NuTo::IntegrationType2D4NConstVariableIp::GetIntegrationPointWeight(int) const
{
    return 4./mNumIp;
}


std::string NuTo::IntegrationType2D4NConstVariableIp::GetStrIdentifier()const
{
    return "2D4NConst"+ std::to_string(mNumIp) + "Ip";
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
