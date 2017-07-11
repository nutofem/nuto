#include "mechanics/integrationtypes/IntegrationTypeEnum.h"
#include <map>
#include <algorithm>
#include "mechanics/MechanicsException.h"

namespace NuTo
{

std::map<eIntegrationType, std::string> GetIntegrationTypeMap()
{
    std::map<eIntegrationType, std::string> map = {
            {eIntegrationType::IntegrationType0DBoundary, "INTEGRATIONTYPE0DBOUNDARY"},
            {eIntegrationType::IntegrationType1D2NGauss1Ip, "INTEGRATIONTYPE1D2NGAUSS1IP"},
            {eIntegrationType::IntegrationType1D2NGauss2Ip, "INTEGRATIONTYPE1D2NGAUSS2IP"},
            {eIntegrationType::IntegrationType1D2NBoundaryGauss3Ip, "INTEGRATIONTYPE1D2NBOUNDARYGAUSS3IP"},
            {eIntegrationType::IntegrationType1D2NGauss3Ip, "INTEGRATIONTYPE1D2NGAUSS3IP"},
            {eIntegrationType::IntegrationType1D2NGauss4Ip, "INTEGRATIONTYPE1D2NGAUSS4IP"},
            {eIntegrationType::IntegrationType1D2NGauss5Ip, "INTEGRATIONTYPE1D2NGAUSS5IP"},
            {eIntegrationType::IntegrationType1D2NLobatto3Ip, "INTEGRATIONTYPE1D2NLOBATTO3IP"},
            {eIntegrationType::IntegrationType1D2NLobatto4Ip, "INTEGRATIONTYPE1D2NLOBATTO4IP"},
            {eIntegrationType::IntegrationType1D2NLobatto5Ip, "INTEGRATIONTYPE1D2NLOBATTO5IP"},
            {eIntegrationType::IntegrationType1D2NLobatto6Ip, "INTEGRATIONTYPE1D2NLOBATTO6IP"},
            {eIntegrationType::IntegrationType2D3NGauss1Ip, "INTEGRATIONTYPE2D3NGAUSS1IP"},
            {eIntegrationType::IntegrationType2D3NGauss3Ip, "INTEGRATIONTYPE2D3NGAUSS3IP"},
            {eIntegrationType::IntegrationType2D3NGauss4Ip, "INTEGRATIONTYPE2D3NGAUSS4IP"},
            {eIntegrationType::IntegrationType2D3NGauss6Ip, "INTEGRATIONTYPE2D3NGAUSS6IP"},
            {eIntegrationType::IntegrationType2D3NGauss12Ip, "INTEGRATIONTYPE2D3NGAUSS12IP"},
            {eIntegrationType::IntegrationType2D3NGauss12IpDetail, "INTEGRATIONTYPE2D3NGAUSS12IPDETAIL"},
            {eIntegrationType::IntegrationType2D3NGauss13Ip, "INTEGRATIONTYPE2D3NGAUSS13IP"},
            {eIntegrationType::IntegrationType2D3NGauss16Ip, "INTEGRATIONTYPE2D3NGAUSS16IP"},
            {eIntegrationType::IntegrationType2D4NGauss1Ip, "INTEGRATIONTYPE2D4NGAUSS1IP"},
            {eIntegrationType::IntegrationType2D4NGauss4Ip, "INTEGRATIONTYPE2D4NGAUSS4IP"},
            {eIntegrationType::IntegrationType2D4NGauss9Ip, "INTEGRATIONTYPE2D4NGAUSS9IP"},
            {eIntegrationType::IntegrationType2D4NLobatto9Ip, "INTEGRATIONTYPE2D4NLOBATTO9IP"},
            {eIntegrationType::IntegrationType2D4NLobatto16Ip, "INTEGRATIONTYPE2D4NLOBATTO16IP"},
            {eIntegrationType::IntegrationType2D4NLobatto25Ip, "INTEGRATIONTYPE2D4NLOBATTO25IP"},
            {eIntegrationType::IntegrationType3D4NGauss1Ip, "INTEGRATIONTYPE3D4NGAUSS1IP"},
            {eIntegrationType::IntegrationType3D4NGauss4Ip, "INTEGRATIONTYPE3D4NGAUSS4IP"},
            {eIntegrationType::IntegrationType3D6NGauss1Ip, "INTEGRATIONTYPE3D6NGAUSS1IP"},
            {eIntegrationType::IntegrationType3D6NGauss2x3Ip, "INTEGRATIONTYPE3D6NGAUSS2X3IP"},
            {eIntegrationType::IntegrationType3D8NGauss1Ip, "INTEGRATIONTYPE3D8NGAUSS1IP"},
            {eIntegrationType::IntegrationType3D8NGauss2x2x2Ip, "INTEGRATIONTYPE3D8NGAUSS2X2X2IP"},
            {eIntegrationType::IntegrationType3D8NLobatto3x3x3Ip, "INTEGRATIONTYPE3D8NLOBATTO3X3X3IP"},
            {eIntegrationType::IntegrationType3D8NLobatto4x4x4Ip, "INTEGRATIONTYPE3D8NLOBATTO4X4X4IP"},
            {eIntegrationType::IntegrationType3D8NLobatto5x5x5Ip, "INTEGRATIONTYPE3D8NLOBATTO5X5X5IP"},
            {eIntegrationType::NotSet, ""}};
    return map;
}

std::string IntegrationTypeToString(const eIntegrationType& rIntegrationType)
{
    try
    {
        return GetIntegrationTypeMap().find(rIntegrationType)->second;
    }
    catch (const std::out_of_range& e)
    {
        throw MechanicsException(__PRETTY_FUNCTION__, "Enum undefined or not implemented.");
    }
}

eIntegrationType IntegrationTypeToEnum(std::string rIntegrationType)
{
    std::transform(rIntegrationType.begin(), rIntegrationType.end(), rIntegrationType.begin(), ::toupper);

    for (auto& entry : GetIntegrationTypeMap())
        if (entry.second == rIntegrationType)
            return entry.first;

    throw MechanicsException(__PRETTY_FUNCTION__,
                             "IntegrationType " + rIntegrationType + " has no enum equivalent or is not implemented.");
}

} // namespace NuTo
