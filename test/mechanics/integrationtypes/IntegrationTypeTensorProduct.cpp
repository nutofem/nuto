#include "BoostUnitTest.h"
#include "mechanics/integrationtypes/IntegrationTypeTensorProduct.h"
#include <vector>
#include <iostream>

BOOST_AUTO_TEST_CASE(NodesAndWeightsLobatto)
{
    std::vector<std::vector<double>> nodesReference = {
            {0.},
            {-sqrt(1. / 3.), sqrt(1. / 3.)},
            {-sqrt(3. / 5.), 0., sqrt(3. / 5.)},
            {-sqrt(3. / 7. + 2. / 7. * sqrt(6. / 5.)), -sqrt(3. / 7. - 2. / 7. * sqrt(6. / 5.)),
             sqrt(3. / 7. - 2. / 7. * sqrt(6. / 5.)), sqrt(3. / 7. + 2. / 7. * sqrt(6. / 5.))}};

    std::vector<std::vector<double>> weightsReference = {
            {2},
            {1., 1.},
            {5. / 9., 8. / 9., 5. / 9.},
            {(18. - sqrt(30)) / 36., (18. + sqrt(30)) / 36., (18. + sqrt(30)) / 36., (18. - sqrt(30)) / 36.},
    };

    for (size_t j = 0; j < nodesReference.size(); j++)
    {
        int numIps = j + 1;
        std::vector<double> nodes(numIps);
        std::vector<double> weights(numIps);

        NuTo::IntegrationTypeTensorProduct<1> intType(numIps, NuTo::eIntegrationMethod::GAUSS);
        for (int i = 0; i < numIps; i++)
        {
            nodes[i] = intType.GetLocalIntegrationPointCoordinates(i)(0, 0);
            weights[i] = intType.GetIntegrationPointWeight(i);
        }

        for (size_t i = 0; i < nodes.size(); i++)
        {
            BOOST_CHECK_CLOSE(nodes[i], nodesReference[j][i], 1.e-10);
            BOOST_CHECK_CLOSE(weights[i], weightsReference[j][i], 1.e-10);
        }
    }
}

BOOST_AUTO_TEST_CASE(NodesAndWeightsGauss)
{
    std::vector<std::vector<double>> nodesReference = {
            {0.},
            {-sqrt(1. / 3.), sqrt(1. / 3.)},
            {-sqrt(3. / 5.), 0., sqrt(3. / 5.)},
            {-sqrt(3. / 7. + 2. / 7. * sqrt(6. / 5.)), -sqrt(3. / 7. - 2. / 7. * sqrt(6. / 5.)),
             sqrt(3. / 7. - 2. / 7. * sqrt(6. / 5.)), sqrt(3. / 7. + 2. / 7. * sqrt(6. / 5.))}};

    std::vector<std::vector<double>> weightsReference = {
            {2},
            {1., 1.},
            {5. / 9., 8. / 9., 5. / 9.},
            {(18. - sqrt(30)) / 36., (18. + sqrt(30)) / 36., (18. + sqrt(30)) / 36., (18. - sqrt(30)) / 36.},
    };

    for (size_t j = 0; j < nodesReference.size(); j++)
    {
        int numIps = j + 1;
        std::vector<double> nodes(numIps);
        std::vector<double> weights(numIps);

        NuTo::IntegrationTypeTensorProduct<1> intType(numIps, NuTo::eIntegrationMethod::GAUSS);
        for (int i = 0; i < numIps; i++)
        {
            nodes[i] = intType.GetLocalIntegrationPointCoordinates(i)(0, 0);
            weights[i] = intType.GetIntegrationPointWeight(i);
        }

        for (size_t i = 0; i < nodes.size(); i++)
        {
            BOOST_CHECK_CLOSE(nodes[i], nodesReference[j][i], 1.e-10);
            BOOST_CHECK_CLOSE(weights[i], weightsReference[j][i], 1.e-10);
        }
    }
}

double integrate(std::function<double(Eigen::VectorXd)> f, NuTo::IntegrationTypeBase &intType) {
    double result = 0.;
    for (int i = 0; i < intType.GetNumIntegrationPoints(); i++) {
        double y = f(intType.GetLocalIntegrationPointCoordinates(i));
        double w = intType.GetIntegrationPointWeight(i);
        result += w * y;
    }
    return (result);
}

void Integrate1DLine(int maxOrder, int numNodes1D, NuTo::eIntegrationMethod method) {
    NuTo::IntegrationTypeTensorProduct<1> intType(numNodes1D, method);
    for (int i=0; i<=maxOrder; i++) {
        auto f = [i](Eigen::VectorXd x){ return( std::pow(x[0],i) ); };
        double computedResult = integrate(f, intType);
        double expectedResult = 1./(i+1) * (1. -  std::pow(-1,i+1));
        BOOST_CHECK_SMALL(computedResult - expectedResult, 1.e-13);
    }
}

void Integrate2DQuad(int maxOrder, int numNodes1D, NuTo::eIntegrationMethod method) {
    NuTo::IntegrationTypeTensorProduct<2> intType(numNodes1D, method);
    for (int n=0; n<=maxOrder; n++) {
        for (int i = 0; i < n + 1; i++) {
            auto f = [n,i](Eigen::VectorXd x){ return( std::pow(x[0],i) * std::pow(x[1],n-i) ); };
            double computedResult = integrate(f, intType);
            double expectedResult = 1./(i+1)/(n-i+1) * (   (1. - std::pow(-1,i+1) )*(1. - std::pow(-1,n-i+1) ) );
            BOOST_CHECK_SMALL(computedResult - expectedResult, 1.e-13);
        }
    }
}

void Integrate3DBrick(int maxOrder, int numNodes1D, NuTo::eIntegrationMethod method) {
    NuTo::IntegrationTypeTensorProduct<3> intType(numNodes1D, method);
    for (int n=0; n<=maxOrder; n++) {
        for (int i = 0; i < n + 1; i++) {
            for (int j = 0; j < (n + 1 - i); j++) {
                auto f = [n,i,j](Eigen::VectorXd x){ return( std::pow(x[0],i) * std::pow(x[1],j) * std::pow(x[2],n-i-j) ); };
                double computedResult = integrate(f, intType);
                double expectedResult = 1./(i+1)/(j+1)/(n-i-j+1) * (   (1. - std::pow(-1,i+1) )*(1. - std::pow(-1,j+1) ) * (1. - std::pow(-1,n-i-j+1)) );
                BOOST_CHECK_SMALL(computedResult - expectedResult, 1.e-13);
            }
        }
    }
}

BOOST_AUTO_TEST_CASE(Integrate1DLineGauss)
{
    for (int numNodes1D = 1; numNodes1D < 8; numNodes1D ++) {
        // Gauss integration should be exact for polynomials up to order 2n-1
        Integrate1DLine(2*numNodes1D-1, numNodes1D, NuTo::eIntegrationMethod::GAUSS);
    }
}

BOOST_AUTO_TEST_CASE(Integrate2DQuadGauss)
{
    for (int numNodes1D = 1; numNodes1D < 8; numNodes1D ++) {
        // Gauss integration should be exact for polynomials up to order 2n-1
        Integrate2DQuad(2*numNodes1D-1, numNodes1D, NuTo::eIntegrationMethod::GAUSS);
    }
}

BOOST_AUTO_TEST_CASE(Integrate3DBrickGauss)
{
    for (int numNodes1D = 1; numNodes1D < 8; numNodes1D ++) {
        // Gauss integration should be exact for polynomials up to order 2n-1
        Integrate3DBrick(2*numNodes1D-1, numNodes1D, NuTo::eIntegrationMethod::GAUSS);
    }
}

BOOST_AUTO_TEST_CASE(Integrate1DLineLobatto)
{
    for (int numNodes1D = 2; numNodes1D < 8; numNodes1D ++) {
        // Lobatto integration should be exact for polynomials up to order 2n-3
        Integrate1DLine(2*numNodes1D-3, numNodes1D, NuTo::eIntegrationMethod::LOBATTO);
    }
}

BOOST_AUTO_TEST_CASE(Integrate2DQuadLobatto)
{
    for (int numNodes1D = 2; numNodes1D < 8; numNodes1D ++) {
        // Lobatto integration should be exact for polynomials up to order 2n-3
        Integrate2DQuad(2*numNodes1D-3, numNodes1D, NuTo::eIntegrationMethod::LOBATTO);
    }
}

BOOST_AUTO_TEST_CASE(Integrate3DBrickLobatto)
{
    for (int numNodes1D = 2; numNodes1D < 8; numNodes1D ++) {
        // Lobatto integration should be exact for polynomials up to order 2n-3
        Integrate3DBrick(2*numNodes1D-3, numNodes1D, NuTo::eIntegrationMethod::LOBATTO);
    }
}
