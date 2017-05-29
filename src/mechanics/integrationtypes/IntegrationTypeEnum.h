// $Id$ 
#pragma once

#include <string>

namespace NuTo
{

//! @author Joerg F. Unger
//! @date May 2, 2010
//! @brief ...

enum class eIntegrationType
{
    IntegrationType0DBoundary,
    IntegrationType1D2NGauss1Ip,
    IntegrationType1D2NGauss2Ip,
    IntegrationType1D2NBoundaryGauss3Ip,
    IntegrationType1D2NGauss3Ip,
    IntegrationType1D2NGauss4Ip,
    IntegrationType1D2NGauss5Ip,
    IntegrationType1D2NLobatto3Ip,
    IntegrationType1D2NLobatto4Ip,
    IntegrationType1D2NLobatto5Ip,
    IntegrationType1D2NLobatto6Ip,
    IntegrationType2D3NGauss1Ip,
    IntegrationType2D3NGauss3Ip,
    IntegrationType2D3NGauss4Ip,
    IntegrationType2D3NGauss6Ip,
    IntegrationType2D3NGauss12Ip,
    IntegrationType2D3NGauss12IpDetail,
    IntegrationType2D3NGauss13Ip,
    IntegrationType2D3NGauss16Ip,
    IntegrationType2D4NGauss1Ip,
    IntegrationType2D4NGauss4Ip,
    IntegrationType2D4NGauss9Ip,
    IntegrationType2D4NLobatto9Ip,
    IntegrationType2D4NLobatto16Ip,
    IntegrationType2D4NLobatto25Ip,
    IntegrationType3D4NGauss1Ip,
    IntegrationType3D4NGauss4Ip,
    IntegrationType3D6NGauss1Ip,
    IntegrationType3D6NGauss2x3Ip,
    IntegrationType3D8NGauss1Ip,
    IntegrationType3D8NGauss2x2x2Ip,
    IntegrationType3D8NLobatto3x3x3Ip,
    IntegrationType3D8NLobatto4x4x4Ip,
    IntegrationType3D8NLobatto5x5x5Ip,
    NotSet
};

std::string IntegrationTypeToString(const eIntegrationType& rIntegrationType);
eIntegrationType IntegrationTypeToEnum(std::string rIntegrationType);

}// namespace NuTo
