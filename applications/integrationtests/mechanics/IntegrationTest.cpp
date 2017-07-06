#include "BoostUnitTest.h"

#include "mechanics/integrationtypes/IntegrationType2D.h"
#include "mechanics/integrationtypes/IntegrationType2D4NGauss1Ip.h"
#include "mechanics/integrationtypes/IntegrationType2D4NGauss4Ip.h"
#include "mechanics/integrationtypes/IntegrationType2D4NGauss9Ip.h"
#include "mechanics/integrationtypes/IntegrationType2D4NLobatto16Ip.h"
#include "mechanics/integrationtypes/IntegrationType2D4NLobatto25Ip.h"
#include "mechanics/integrationtypes/IntegrationType2D4NLobatto9Ip.h"

#include "mechanics/integrationtypes/IntegrationType2D3NGauss1Ip.h"
#include "mechanics/integrationtypes/IntegrationType2D3NGauss3Ip.h"
#include "mechanics/integrationtypes/IntegrationType2D3NGauss4Ip.h"
#include "mechanics/integrationtypes/IntegrationType2D3NGauss6Ip.h"
#include "mechanics/integrationtypes/IntegrationType2D3NGauss12Ip.h"
#include "mechanics/integrationtypes/IntegrationType2D3NGauss12IpDetail.h"
#include "mechanics/integrationtypes/IntegrationType2D3NGauss13Ip.h"
#include "mechanics/integrationtypes/IntegrationType2D3NGauss16Ip.h"

#include "mechanics/integrationtypes/IntegrationType3D4NGauss1Ip.h"
#include "mechanics/integrationtypes/IntegrationType3D4NGauss4Ip.h"

#include "mechanics/integrationtypes/IntegrationType3D8NGauss1Ip.h"
#include "mechanics/integrationtypes/IntegrationType3D8NGauss2x2x2Ip.h"
#include "mechanics/integrationtypes/IntegrationType3D8NLobatto.h"
#include "mechanics/integrationtypes/IntegrationType3D8NLobatto_Def.h"

#include <boost/filesystem.hpp>
#include <iostream>

double testPoly2D(Eigen::Vector2d coord, int order) {
    double result = 0;
    double x = coord(0);
    double y = coord(1);
    for (int n = 0; n < order + 1; n++) {
        for (int i = 0; i < n + 1; i++) {
            result += std::pow(x, i) * std::pow(y, n - i);
        }
    }
    return (result);
}

double testPoly3D(Eigen::Vector3d coord, int order) {
    double result = 0;
    double x = coord(0);
    double y = coord(1);
    double z = coord(2);
    for (int n = 0; n < order + 1; n++) {
        for (int i = 0; i < n + 1; i++) {
            for (int j = 0; j < (n + 1 - i); j++) {
                result += std::pow(x, i) * std::pow(y, j) * std::pow(z, n - i - j);
            }
        }
    }
    return (result);
}

double integrate(std::function<double(Eigen::VectorXd)> f, NuTo::IntegrationTypeBase &intType) {
    double result = 0.;
    for (int i = 0; i < intType.GetNumIntegrationPoints(); i++) {
        double y = f(intType.GetLocalIntegrationPointCoordinates(i));
        result += intType.GetIntegrationPointWeight(i) * y;
    }
    return (result);
}

double integrate2D(int order, NuTo::IntegrationType2D &intType) {
    std::function<double(Eigen::VectorXd)> f = [order](Eigen::VectorXd x) {
        return(testPoly2D(x,order));
    };
    return(integrate(f,intType));
}

double integrate3D(int order, NuTo::IntegrationType3D &intType) {
    std::function<double(Eigen::VectorXd)> f = [order](Eigen::VectorXd x) {
        return(testPoly3D(x,order));
    };
    return(integrate(f,intType));
}

BOOST_AUTO_TEST_CASE(PolynomialIntegrationQuad) {
    // 2D
    NuTo::IntegrationType2D4NGauss1Ip intGaussQuad1;
    NuTo::IntegrationType2D4NGauss4Ip intGaussQuad2;
    NuTo::IntegrationType2D4NGauss9Ip intGaussQuad3;

    NuTo::IntegrationType2D4NLobatto9Ip intLobQuad3;
    NuTo::IntegrationType2D4NLobatto16Ip intLobQuad4;
    NuTo::IntegrationType2D4NLobatto25Ip intLobQuad5;

    std::vector<NuTo::IntegrationType2D *> intTypesQuad2D = {
        &intGaussQuad1, &intGaussQuad2, &intGaussQuad3,
        &intLobQuad3,   &intLobQuad4,   &intLobQuad5};

    // Analytical results (result[i] = integral of polynomial of order i with all
    // coefficients set to 1)
    std::vector<double> analyticalResultsQuad = {4.,
                                               4.,
                                               20. / 3.,
                                               20. / 3.,
                                               392. / 45,
                                               392. / 45,
                                               3272. / 315,
                                               3272. / 315,
                                               2068. / 175,
                                               2068. / 175};

    int polyOrder = 5;

    std::cout << "----------------------------" << std::endl;
    std::cout << "Expected Result Quad: " << analyticalResultsQuad[polyOrder]
                 << std::endl;

    for (NuTo::IntegrationType2D *intType : intTypesQuad2D) {
        std::cout << "Int Points " << intType->GetNumIntegrationPoints()
                  << "Result: " << integrate2D(polyOrder, *(intType)) << std::endl;
    }
}

BOOST_AUTO_TEST_CASE(PolynomialIntegrationTriangle) {
    NuTo::IntegrationType2D3NGauss1Ip intGaussTri1;
    NuTo::IntegrationType2D3NGauss3Ip intGaussTri3;
    NuTo::IntegrationType2D3NGauss4Ip intGaussTri4;
    NuTo::IntegrationType2D3NGauss6Ip intGaussTri6;
    NuTo::IntegrationType2D3NGauss12Ip intGaussTri12;
    NuTo::IntegrationType2D3NGauss12IpDetail intGaussTri12d;
    NuTo::IntegrationType2D3NGauss13Ip intGaussTri13;
    NuTo::IntegrationType2D3NGauss16Ip intGaussTri16;

    std::vector<NuTo::IntegrationType2D *> intTypesTri2D = {
        &intGaussTri1, &intGaussTri3, &intGaussTri4, &intGaussTri6,
        &intGaussTri12, &intGaussTri12d, &intGaussTri13, &intGaussTri16
    };

    // Analytical results (result[i] = integral of polynomial of order i with all
    // coefficients set to 1)
    std::vector<double> analyticalResultsTriangle = {
        1./2.,
        5./6.,
        25./24.,
        47./40.,
        91./72.,
        3341./2520.,
        13817./10080.
    };

    int polyOrder = 5;

    std::cout << "----------------------------" << std::endl;
    std::cout << "Expected Result Triangle: " << analyticalResultsTriangle[polyOrder]
                 << std::endl;

    for (NuTo::IntegrationType2D *intType : intTypesTri2D) {
        std::cout << "Int Points " << intType->GetNumIntegrationPoints()
                  << "Result: " << integrate2D(polyOrder, *(intType)) << std::endl;
    }
}

BOOST_AUTO_TEST_CASE(PolynomialIntegrationTet) {
    // Tetrahedron
    NuTo::IntegrationType3D4NGauss1Ip intGaussTet1;
    NuTo::IntegrationType3D4NGauss4Ip intGaussTet4;

    std::vector<NuTo::IntegrationType3D *> intTypesTet = {
        &intGaussTet1, &intGaussTet4
    };

    // Analytical results (result[i] = integral of polynomial of order i with all
    // coefficients set to 1)
    std::vector<double> analyticalResultsTet = {
            1./6.,
            7./24.,
            11./30.,
            59./144.,
            313./720.,
            9067./20160.,
            83317181440.
    };

    int polyOrder = 5;

    std::cout << "----------------------------" << std::endl;
    std::cout << "Expected Result Tet: " << analyticalResultsTet[polyOrder]
                 << std::endl;

    for (NuTo::IntegrationType3D *intType : intTypesTet) {
        std::cout << "Int Points " << intType->GetNumIntegrationPoints()
                  << "Result: " << integrate3D(polyOrder, *(intType)) << std::endl;
    }
}

BOOST_AUTO_TEST_CASE(PolynomialIntegrationBrick) {
    // 3D
    NuTo::IntegrationType3D8NGauss1Ip intGaussBrick1;
    NuTo::IntegrationType3D8NGauss2x2x2Ip intGaussBrick2;

    NuTo::IntegrationType3D8NLobatto<3> intLobBrick3;
    NuTo::IntegrationType3D8NLobatto<4> intLobBrick4;
    NuTo::IntegrationType3D8NLobatto<5> intLobBrick5;

    std::vector<NuTo::IntegrationType3D *> intTypesBrick3D = {
        &intGaussBrick1, &intGaussBrick2, &intLobBrick3, &intLobBrick4, &intLobBrick5};

    // Analytical results (result[i] = integral of polynomial of order i with all
    // coefficients set to 1)
    std::vector<double> analyticalResultsBrick = {8.,
                                               8.,
                                               16.,
                                               16.,
                                               352. / 15,
                                               352. / 15,
                                               5744. / 189,
                                               5744. / 189,
                                               174056. / 4725,
                                               174056. / 4725};

    int polyOrder = 5;

    std::cout << "----------------------------" << std::endl;
    std::cout << "Expected Result Brick: " << analyticalResultsBrick[polyOrder]
                 << std::endl;

    for (NuTo::IntegrationType3D *intType : intTypesBrick3D) {
        std::cout << "Int Points " << intType->GetNumIntegrationPoints()
                  << "Result: " << integrate3D(polyOrder, *(intType)) << std::endl;
    }
}
