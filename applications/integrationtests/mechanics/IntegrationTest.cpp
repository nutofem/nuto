#include "BoostUnitTest.h"

#include "mechanics/integrationtypes/IntegrationType2D.h"
#include "mechanics/integrationtypes/IntegrationType2D4NGauss1Ip.h"
#include "mechanics/integrationtypes/IntegrationType2D4NGauss4Ip.h"
#include "mechanics/integrationtypes/IntegrationType2D4NGauss9Ip.h"
#include "mechanics/integrationtypes/IntegrationType2D4NLobatto16Ip.h"
#include "mechanics/integrationtypes/IntegrationType2D4NLobatto25Ip.h"
#include "mechanics/integrationtypes/IntegrationType2D4NLobatto9Ip.h"

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

double integrate2D(int order, NuTo::IntegrationType2D &intType) {
  int N = intType.GetNumIntegrationPoints();
  std::vector<double> weights(N);
  std::vector<Eigen::VectorXd> points(N);
  for (int i = 0; i < N; i++) {
    weights[i] = intType.GetIntegrationPointWeight(i);
    points[i] = intType.GetLocalIntegrationPointCoordinates(i);
  }
  double result = 0.;
  for (int i = 0; i < N; i++) {
    result += weights[i] * testPoly2D(points[i], order);
  }
  return (result);
}

double integrate3D(int order, NuTo::IntegrationType3D &intType) {
  int N = intType.GetNumIntegrationPoints();
  std::vector<double> weights(N);
  std::vector<Eigen::VectorXd> points(N);
  for (int i = 0; i < N; i++) {
    weights[i] = intType.GetIntegrationPointWeight(i);
    points[i] = intType.GetLocalIntegrationPointCoordinates(i);
  }
  double result = 0.;
  for (int i = 0; i < N; i++) {
    result += weights[i] * testPoly3D(points[i], order);
  }
  return (result);
}

BOOST_AUTO_TEST_CASE(PolynomialIntegration2D) {
  // 2D
  NuTo::IntegrationType2D4NGauss1Ip *intGaussQuad1 =
      new NuTo::IntegrationType2D4NGauss1Ip();
  NuTo::IntegrationType2D4NGauss4Ip *intGaussQuad2 =
      new NuTo::IntegrationType2D4NGauss4Ip();
  NuTo::IntegrationType2D4NGauss9Ip *intGaussQuad3 =
      new NuTo::IntegrationType2D4NGauss9Ip();

  NuTo::IntegrationType2D4NLobatto9Ip *intLobQuad3 =
      new NuTo::IntegrationType2D4NLobatto9Ip();
  NuTo::IntegrationType2D4NLobatto16Ip *intLobQuad4 =
      new NuTo::IntegrationType2D4NLobatto16Ip();
  NuTo::IntegrationType2D4NLobatto25Ip *intLobQuad5 =
      new NuTo::IntegrationType2D4NLobatto25Ip();

  std::vector<NuTo::IntegrationType2D *> intTypesQuad2D = {
      intGaussQuad1, intGaussQuad2, intGaussQuad3,
      intLobQuad3,   intLobQuad4,   intLobQuad5};

  // Analytical results (result[i] = integral of polynomial of order i with all
  // coefficients set to 1)
  std::vector<double> analyticalResults2D = {4.,
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

  std::cout << "Expected Result 2D: " << analyticalResults2D[polyOrder]
            << std::endl;

  for (NuTo::IntegrationType2D *intType : intTypesQuad2D) {
    std::cout << "Int Points " << intType->GetNumIntegrationPoints()
              << "Result: " << integrate2D(polyOrder, *(intType)) << std::endl;
  }

  delete intGaussQuad1;
  delete intGaussQuad2;
  delete intGaussQuad3;
  delete intLobQuad3;
  delete intLobQuad4;
  delete intLobQuad5;
}

BOOST_AUTO_TEST_CASE(PolynomialIntegration3D) {
  // 3D
  NuTo::IntegrationType3D8NGauss1Ip *intGaussBrick1 =
      new NuTo::IntegrationType3D8NGauss1Ip();
  NuTo::IntegrationType3D8NGauss2x2x2Ip *intGaussBrick2 =
      new NuTo::IntegrationType3D8NGauss2x2x2Ip();

  NuTo::IntegrationType3D8NLobatto<3> *intLobBrick3 =
      new NuTo::IntegrationType3D8NLobatto<3>();
  NuTo::IntegrationType3D8NLobatto<4> *intLobBrick4 =
      new NuTo::IntegrationType3D8NLobatto<4>();
  NuTo::IntegrationType3D8NLobatto<5> *intLobBrick5 =
      new NuTo::IntegrationType3D8NLobatto<5>();

  std::vector<NuTo::IntegrationType3D *> intTypesBrick3D = {
      intGaussBrick1, intGaussBrick2, intLobBrick3, intLobBrick4, intLobBrick5};

  // Analytical results (result[i] = integral of polynomial of order i with all
  // coefficients set to 1)
  std::vector<double> analyticalResults3D = {8.,
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

  std::cout << "Expected Result 3D: " << analyticalResults3D[polyOrder]
            << std::endl;

  for (NuTo::IntegrationType3D *intType : intTypesBrick3D) {
    std::cout << "Int Points " << intType->GetNumIntegrationPoints()
              << "Result: " << integrate3D(polyOrder, *(intType)) << std::endl;
  }

  delete intGaussBrick1;
  delete intGaussBrick2;
  delete intLobBrick3;
  delete intLobBrick4;
  delete intLobBrick5;
}
