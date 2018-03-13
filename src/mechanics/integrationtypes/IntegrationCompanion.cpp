#include "IntegrationCompanion.h"

#include "base/Exception.h"

#include "IntegrationType2D3NGauss1Ip.h"
#include "IntegrationType2D3NGauss3Ip.h"
#include "IntegrationType2D3NGauss4Ip.h"
#include "IntegrationType2D3NGauss6Ip.h"
#include "IntegrationType2D3NGauss12Ip.h"
#include "IntegrationType2D3NGauss13Ip.h"
#include "IntegrationType2D3NGauss16Ip.h"
#include "IntegrationTypeTriangle.h"

#include "IntegrationType3D4NGauss1Ip.h"
#include "IntegrationType3D4NGauss4Ip.h"
#include "IntegrationTypeTetrahedron.h"

#include "IntegrationType3D6NGauss1Ip.h"
#include "IntegrationType3D6NGauss2x3Ip.h"

#include "IntegrationTypeTensorProduct.h"

using namespace NuTo;

namespace
{

std::unique_ptr<IntegrationTypeBase> GaussLine(int order)
{
    return std::make_unique<IntegrationTypeTensorProduct<1>>(order, eIntegrationMethod::GAUSS);
}


std::unique_ptr<IntegrationTypeBase> GaussTriangle(int order)
{
    switch (order)
    {
    case 1:
        return std::make_unique<IntegrationType2D3NGauss1Ip>();
    case 2:
        return std::make_unique<IntegrationType2D3NGauss3Ip>();
    case 3:
        return std::make_unique<IntegrationType2D3NGauss4Ip>();
    case 4:
        return std::make_unique<IntegrationType2D3NGauss6Ip>();
    case 5:
    case 6:
        return std::make_unique<IntegrationType2D3NGauss12Ip>();
    case 7:
    // return std::make_unique<IntegrationType2D3NGauss13Ip>();
    case 8:
    case 9:
    case 10:
    case 11:
    case 12:
    case 13:
    case 14:
    case 15:
    case 16:
    case 17:
    case 18:
    case 19:
    case 20:
        return std::make_unique<IntegrationTypeTriangle>(order);
    default:
        throw Exception(__PRETTY_FUNCTION__,
                        "Triangle integration of order " + std::to_string(order) + " is not defined.");
    }
}


std::unique_ptr<IntegrationTypeBase> GaussQuadrilateral(int order)
{
    return std::make_unique<IntegrationTypeTensorProduct<2>>(order, eIntegrationMethod::GAUSS);
}


std::unique_ptr<IntegrationTypeBase> GaussTetrahedron(int order)
{
    switch (order)
    {
    case 1:
        return std::make_unique<IntegrationType3D4NGauss1Ip>();
    case 2:
        return std::make_unique<IntegrationType3D4NGauss4Ip>();
    case 3:
    case 4:
    case 5:
    case 6:
    case 7:
    case 8:
    case 9:
    case 10:
    case 11:
    case 12:
    case 13:
    case 14:
    case 15:
    case 16:
    case 17:
    case 18:
    case 19:
    case 20:
        return std::make_unique<IntegrationTypeTetrahedron>(order);
    default:
        throw Exception(__PRETTY_FUNCTION__, "Tet integration of order " + std::to_string(order) + " is not defined.");
    }
}


std::unique_ptr<IntegrationTypeBase> GaussHex(int order)
{
    return std::make_unique<IntegrationTypeTensorProduct<3>>(order, eIntegrationMethod::GAUSS);
}


std::unique_ptr<IntegrationTypeBase> GaussPrism(int order)
{
    switch (order)
    {
    case 1:
        return std::make_unique<IntegrationType3D6NGauss1Ip>();
    case 2:
        return std::make_unique<IntegrationType3D6NGauss2x3Ip>();
    default:
        throw Exception(__PRETTY_FUNCTION__,
                        "Prism integration of order " + std::to_string(order) + " is not defined.");
    }
}


std::unique_ptr<IntegrationTypeBase> GaussPyramid(int order)
{
    switch (order)
    {
    default:
        throw Exception(__PRETTY_FUNCTION__,
                        "Prism integration of order " + std::to_string(order) + " is not defined.");
    }
}

std::unique_ptr<IntegrationTypeBase> LobattoLine(int order)
{
    return std::make_unique<IntegrationTypeTensorProduct<1>>(order, eIntegrationMethod::LOBATTO);
}

std::unique_ptr<IntegrationTypeBase> LobattoQuadrilateral(int order)
{
    return std::make_unique<IntegrationTypeTensorProduct<2>>(order, eIntegrationMethod::LOBATTO);
}

std::unique_ptr<IntegrationTypeBase> LobattoHex(int order)
{
    return std::make_unique<IntegrationTypeTensorProduct<3>>(order, eIntegrationMethod::LOBATTO);
}
}

std::unique_ptr<IntegrationTypeBase> NuTo::CreateGaussIntegrationType(const Shape& shape, int order)
{
    switch (shape.Enum())
    {
    case eShape::Line:
        return GaussLine(order);
    case eShape::Triangle:
        return GaussTriangle(order);
    case eShape::Quadrilateral:
        return GaussQuadrilateral(order);
    case eShape::Tetrahedron:
        return GaussTetrahedron(order);
    case eShape::Hexahedron:
        return GaussHex(order);
    case eShape::Prism:
        return GaussPrism(order);
    case eShape::Pyramid:
        return GaussPyramid(order);
    default:
        throw Exception(__PRETTY_FUNCTION__, "Shape enum is not known.");
    }
}

std::unique_ptr<IntegrationTypeBase> NuTo::CreateLobattoIntegrationType(const Shape& shape, int order)
{
    switch (shape.Enum())
    {
    case eShape::Line:
        return LobattoLine(order);
    case eShape::Quadrilateral:
        return LobattoQuadrilateral(order);
    case eShape::Hexahedron:
        return LobattoHex(order);
    case eShape::Triangle:
    case eShape::Tetrahedron:
    case eShape::Prism:
    case eShape::Pyramid:
        throw Exception(__PRETTY_FUNCTION__, "Lobatto integration for this shape is not defined.");
    default:
        throw Exception(__PRETTY_FUNCTION__, "Shape enum is not known.");
    }
}
