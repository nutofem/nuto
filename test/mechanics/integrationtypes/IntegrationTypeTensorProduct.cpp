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


