// $Id$
#include <sstream>
#include <iostream>
#include "nuto/mechanics/integrationtypes/IntegrationType1D2NConstVariableIp.h"
#include "nuto/mechanics/MechanicsException.h"
#include <assert.h>

//! @brief constructor
NuTo::IntegrationType1D2NConstVariableIp::IntegrationType1D2NConstVariableIp(int rNumIp)
{
    if (rNumIp<1)
        throw MechanicsException("[NuTo::IntegrationType1D2NConstVariableIp] Number of integration points must be positive.");
    mNumIp = rNumIp;
}

//! @brief returns the local coordinates of an integration point
//! @param rIpNum integration point (counting from zero)
//! @param rCoordinates (result)
void NuTo::IntegrationType1D2NConstVariableIp::GetLocalIntegrationPointCoordinates1D(int rIpNum, double& rCoordinates)const
{
    assert(rIpNum>=0 && rIpNum<mNumIp);
    rCoordinates = -1.+2.*(rIpNum+0.5)/mNumIp;
}

//! @brief returns the total number of integration points for this integration type
//! @return number of integration points
int NuTo::IntegrationType1D2NConstVariableIp::GetNumIntegrationPoints()const
{
    return mNumIp;
}

//! @brief returns the weight of an integration point
//! @param rIpNum integration point (counting from zero)
//! @return weight of integration points
double NuTo::IntegrationType1D2NConstVariableIp::GetIntegrationPointWeight(int rIpNum)const
{
    return 2./mNumIp;
}

//! @brief returns a string with the identifier of the integration type
//! @return identifier
std::string NuTo::IntegrationType1D2NConstVariableIp::GetStrIdentifier()const
{
    std::ostringstream o;
    o << mNumIp;
    return std::string("1D2NConst"+ o.str() + "Ip");
}

#ifdef ENABLE_VISUALIZE
void NuTo::IntegrationType1D2NConstVariableIp::GetVisualizationCells(
    unsigned int& NumVisualizationPoints,
    std::vector<double>& VisualizationPointLocalCoordinates,
    unsigned int& NumVisualizationCells,
    std::vector<NuTo::CellBase::eCellTypes>& VisualizationCellType,
    std::vector<unsigned int>& VisualizationCellsIncidence,
    std::vector<unsigned int>& VisualizationCellsIP) const
{
    throw NuTo::MechanicsException("[NuTo::IntegrationType1D2NConstVariableIp::GetVisualizationCells] Integrationcells are not implemented");
}
#endif // ENABLE_VISUALIZE
