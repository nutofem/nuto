#include "IntegrationCompanion.h"

#include "base/Exception.h"

#include "IntegrationType2D3NGauss1Ip.h"
#include "IntegrationType2D3NGauss3Ip.h"
#include "IntegrationType2D3NGauss4Ip.h"
#include "IntegrationType2D3NGauss6Ip.h"
#include "IntegrationType2D3NGauss12Ip.h"
#include "IntegrationType2D3NGauss13Ip.h"
#include "IntegrationType2D3NGauss16Ip.h"

#include "IntegrationType3D4NGauss1Ip.h"
#include "IntegrationType3D4NGauss4Ip.h"

#include "IntegrationType3D6NGauss1Ip.h"
#include "IntegrationType3D6NGauss2x3Ip.h"

#include "IntegrationTypeTensorProduct.h"

using namespace NuTo;

namespace
{

std::unique_ptr<IntegrationTypeBase> S_Line(int order)
{
    return std::make_unique<IntegrationTypeTensorProduct<1>>(order, eIntegrationMethod::GAUSS);
}


std::unique_ptr<IntegrationTypeBase> S_Triangle(int order)
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
        return std::make_unique<IntegrationType2D3NGauss13Ip>();
    default:
        throw Exception(__PRETTY_FUNCTION__,
                        "Triangle integration of order " + std::to_string(order) + " is not defined.");
    }
}


std::unique_ptr<IntegrationTypeBase> S_Quadrilateral(int order)
{
    return std::make_unique<IntegrationTypeTensorProduct<2>>(order, eIntegrationMethod::GAUSS);
}


std::unique_ptr<IntegrationTypeBase> S_Tetrahedron(int order)
{
    switch (order)
    {
    case 1:
        return std::make_unique<IntegrationType3D4NGauss1Ip>();
    case 2:
        return std::make_unique<IntegrationType3D4NGauss4Ip>();
    default:
        throw Exception(__PRETTY_FUNCTION__, "Tet integration of order " + std::to_string(order) + " is not defined.");
    }
}


std::unique_ptr<IntegrationTypeBase> S_Hex(int order)
{
    return std::make_unique<IntegrationTypeTensorProduct<3>>(order, eIntegrationMethod::GAUSS);
}


std::unique_ptr<IntegrationTypeBase> S_Prism(int order)
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


std::unique_ptr<IntegrationTypeBase> S_Pyramid(int order)
{
    switch (order)
    {
    default:
        throw Exception(__PRETTY_FUNCTION__,
                        "Prism integration of order " + std::to_string(order) + " is not defined.");
    }
}
}

std::unique_ptr<IntegrationTypeBase> NuTo::CreateIntegrationType(const Shape& shape, int order)
{
    switch (shape.Enum())
    {
    case eShape::Line:
        return S_Line(order);
    case eShape::Triangle:
        return S_Triangle(order);
    case eShape::Quadrilateral:
        return S_Quadrilateral(order);
    case eShape::Tetrahedron:
        return S_Tetrahedron(order);
    case eShape::Hexahedron:
        return S_Hex(order);
    case eShape::Prism:
        return S_Prism(order);
    case eShape::Pyramid:
        return S_Pyramid(order);
    default:
        throw Exception(__PRETTY_FUNCTION__, "Shape enum is not known.");
    }
}
