#include "nuto/mechanics/integrationtypes/IntegrationType1D2NConstVariableIp.h"
#include "nuto/mechanics/MechanicsException.h"
#include <assert.h>
#include <string>

NuTo::IntegrationType1D2NConstVariableIp::IntegrationType1D2NConstVariableIp(int rNumIp)
{
    if (rNumIp<1)
        throw MechanicsException(__PRETTY_FUNCTION__, "Number of integration points must be positive.");
    mNumIp = rNumIp;
}


void NuTo::IntegrationType1D2NConstVariableIp::GetLocalIntegrationPointCoordinates1D(int rIpNum, double& rCoordinates) const
{
    assert(rIpNum >= 0 && rIpNum < mNumIp);
    rCoordinates = -2.0*(rIpNum + 0.5)/mNumIp;
}


unsigned int NuTo::IntegrationType1D2NConstVariableIp::GetNumIntegrationPoints() const
{
    return mNumIp;
}


double NuTo::IntegrationType1D2NConstVariableIp::GetIntegrationPointWeight(int) const
{
    return 2./mNumIp;
}


std::string NuTo::IntegrationType1D2NConstVariableIp::GetStrIdentifier()const
{
    return "1D2NConst" + std::to_string(mNumIp) + "Ip";
}

#ifdef ENABLE_VISUALIZE
void NuTo::IntegrationType1D2NConstVariableIp::GetVisualizationCells(
    unsigned int& NumVisualizationPoints,
    std::vector<double>& VisualizationPointLocalCoordinates,
    unsigned int& NumVisualizationCells,
    std::vector<NuTo::eCellTypes>& VisualizationCellType,
    std::vector<unsigned int>& VisualizationCellsIncidence,
    std::vector<unsigned int>& VisualizationCellsIP) const
{
    throw NuTo::MechanicsException("[NuTo::IntegrationType1D2NConstVariableIp::GetVisualizationCells] Integrationcells are not implemented");
}
#endif // ENABLE_VISUALIZE
