// $Id$
#pragma once

namespace NuTo
{

//! @author Joerg F. Unger
//! @date May 2, 2010
//! @brief ...

enum class eIntegrationMethod
{
    GAUSS,
    LOBATTO
};

enum class eIntegrationType
{
    IntegrationType0DBoundary = 0,
    IntegrationType1D2NGauss1Ip,
    IntegrationType1D2NGauss2Ip,
    IntegrationType1D2NBoundaryGauss3Ip,
    IntegrationType1D2NGauss3Ip,
    IntegrationType1D2NGauss4Ip,
    IntegrationType1D2NGauss5Ip,
    IntegrationType1D2NGauss6Ip,
    IntegrationType1D2NGauss7Ip,
    IntegrationType1D2NGauss8Ip,
    IntegrationType1D2NGauss9Ip,
    IntegrationType1D2NGauss10Ip,
    IntegrationType1D2NGauss11Ip,
    IntegrationType1D2NGauss12Ip,
    IntegrationType1D2NGauss13Ip,
    IntegrationType1D2NGauss14Ip,
    IntegrationType1D2NGauss15Ip,
    IntegrationType1D2NGauss16Ip,
    IntegrationType1D2NGauss17Ip,
    IntegrationType1D2NGauss18Ip,
    IntegrationType1D2NGauss19Ip,
    IntegrationType1D2NGauss20Ip,
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
    IntegrationType2D3NLattice3Ip,
    IntegrationType2D4NGauss1Ip,
    IntegrationType2D4NGauss4Ip,
    IntegrationType2D4NGauss9Ip,
    IntegrationTypeTensorProduct2D4N4,
    IntegrationTypeTensorProduct2D4N5,
    IntegrationTypeTensorProduct2D4N6,
    IntegrationTypeTensorProduct2D4N7,
    IntegrationTypeTensorProduct2D4N8,
    IntegrationTypeTensorProduct2D4N9,
    IntegrationTypeTensorProduct2D4N10,
    IntegrationType2D4NLobatto9Ip,
    IntegrationType2D4NLobatto16Ip,
    IntegrationType2D4NLobatto25Ip,
    IntegrationType3D4NGauss1Ip,
    IntegrationType3D4NGauss4Ip,
    IntegrationType3D8NGauss1Ip,
    IntegrationType3D8NGauss2x2x2Ip,
    IntegrationType3D8NLobatto3x3x3Ip,
    IntegrationType3D8NLobatto4x4x4Ip,
    IntegrationType3D8NLobatto5x5x5Ip,
    NumIntegrationTypes,
    NotSet
};

} // namespace NuTo
