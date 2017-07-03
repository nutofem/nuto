#include "BoostUnitTest.h"
#include "mechanics/integrationtypes/IntegrationType1D2NLobatto.h"
#include <vector>
#include <iostream>

BOOST_AUTO_TEST_CASE(NodesAndWeights)
{
    std::vector< std::vector<double> > nodesReference = {
        {},
        {-1.,1.},
        {-1.,0.,1.},
        {-1., -sqrt(5.)/5., sqrt(5.)/5., 1.},
        {-1.,-sqrt(21.)/7.,0.,sqrt(21.)/7.,1.},
        {-1.,-sqrt(1./21.*(7+2.*sqrt(7))), -sqrt(1./21.*(7-2.*sqrt(7))), sqrt(1./21.*(7-2.*sqrt(7))), sqrt(1./21.*(7+2.*sqrt(7))), 1.}
    };

    std::vector< std::vector<double> > weightsReference = {
        {},
        {1.,1.},
        {1./3.,4./3.,1./3.},
        {1./6.,5./6., 5./6., 1./6.},
        {0.1, 49./90., 32./45., 49./90., 0.1},
        {1./15,1./30.*(14.-sqrt(7)),1./30.*(14.+sqrt(7)),1./30.*(14.+sqrt(7)),  1./30.*(14.-sqrt(7)), 1./15}
    };

    int maxOrder = nodesReference.size()-1;
    int minOrder = 2;

    for (int order = minOrder; order<=maxOrder; order++) {
        int numIps = order+1;
        std::vector<double> nodes(numIps);
        std::vector<double> weights(numIps);

        NuTo::IntegrationType1D2NLobatto intType(numIps);
        for (int i=0; i<numIps; i++) {
            nodes[i]   = intType.GetLocalIntegrationPointCoordinates(i)(0,0);
            weights[i] = intType.GetIntegrationPointWeight(i);
        }

        for (size_t i=0; i<nodes.size(); i++) {
            BOOST_CHECK_CLOSE(nodes[i], nodesReference[order][i],1.e-10);
            BOOST_CHECK_CLOSE(weights[i], weightsReference[order][i],1.e-10);
        }
    }

}

