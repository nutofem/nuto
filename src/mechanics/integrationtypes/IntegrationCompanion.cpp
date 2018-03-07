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

std::unique_ptr<IntegrationTypeBase> G_Line(int order)
{
    return std::make_unique<IntegrationTypeTensorProduct<1>>(order, eIntegrationMethod::GAUSS);
}


std::unique_ptr<IntegrationTypeBase> G_Triangle(int order)
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


std::unique_ptr<IntegrationTypeBase> G_Quadrilateral(int order)
{
    return std::make_unique<IntegrationTypeTensorProduct<2>>(order, eIntegrationMethod::GAUSS);
}


std::unique_ptr<IntegrationTypeBase> G_Tetrahedron(int order)
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


std::unique_ptr<IntegrationTypeBase> G_Hex(int order)
{
    return std::make_unique<IntegrationTypeTensorProduct<3>>(order, eIntegrationMethod::GAUSS);
}


std::unique_ptr<IntegrationTypeBase> G_Prism(int order)
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


std::unique_ptr<IntegrationTypeBase> G_Pyramid(int order)
{
    switch (order)
    {
    default:
        throw Exception(__PRETTY_FUNCTION__,
                        "Prism integration of order " + std::to_string(order) + " is not defined.");
    }
}

std::unique_ptr<IntegrationTypeBase> L_Line(int order)
{
    return std::make_unique<IntegrationTypeTensorProduct<1>>(order, eIntegrationMethod::LOBATTO);
}

std::unique_ptr<IntegrationTypeBase> L_Quadrilateral(int order)
{
    return std::make_unique<IntegrationTypeTensorProduct<2>>(order, eIntegrationMethod::LOBATTO);
}

std::unique_ptr<IntegrationTypeBase> L_Hex(int order)
{
    return std::make_unique<IntegrationTypeTensorProduct<3>>(order, eIntegrationMethod::LOBATTO);
}
}

std::unique_ptr<IntegrationTypeBase> NuTo::CreateGaussIntegrationType(const Shape& shape, int order)
{
    switch (shape.Enum())
    {
    case eShape::Line:
        return G_Line(order);
    case eShape::Triangle:
        return G_Triangle(order);
    case eShape::Quadrilateral:
        return G_Quadrilateral(order);
    case eShape::Tetrahedron:
        return G_Tetrahedron(order);
    case eShape::Hexahedron:
        return G_Hex(order);
    case eShape::Prism:
        return G_Prism(order);
    case eShape::Pyramid:
        return G_Pyramid(order);
    default:
        throw Exception(__PRETTY_FUNCTION__, "Shape enum is not known.");
    }
}

std::unique_ptr<IntegrationTypeBase> NuTo::CreateLobattoIntegrationType(const Shape& shape, int order)
{
    switch (shape.Enum())
    {
    case eShape::Line:
        return L_Line(order);
    case eShape::Quadrilateral:
        return L_Quadrilateral(order);
    case eShape::Hexahedron:
        return L_Hex(order);
    case eShape::Triangle:
    case eShape::Tetrahedron:
    case eShape::Prism:
    case eShape::Pyramid:
        throw Exception(__PRETTY_FUNCTION__, "Lobatto integration for this shape is not defined.");
    default:
        throw Exception(__PRETTY_FUNCTION__, "Shape enum is not known.");
    }
}
