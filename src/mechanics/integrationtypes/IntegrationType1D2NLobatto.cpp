#ifdef ENABLE_SERIALIZATION
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/xml_oarchive.hpp>
#include <boost/archive/xml_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#endif //ENABLE_SERIALIZATION

#ifdef ENABLE_VISUALIZE
#include "visualize/VisualizeEnum.h"
#endif // ENABLE_VISUALIZE

#include "mechanics/integrationtypes/IntegrationType1D2NLobatto.h"
#include "mechanics/elements/ElementShapeFunctions.h"
#include "boost/math/special_functions/legendre.hpp"
#include "math/Polynomial.h"

// constructor
NuTo::IntegrationType1D2NLobatto::IntegrationType1D2NLobatto(int numIps)
{
    nIps = numIps;
    std::vector<double> innerNodes = NuTo::Math::Polynomial::legendre_deriv_roots(numIps-1);
    iPts.push_back(-1.);
    for (double value : innerNodes) {
        iPts.push_back(value);
    }
    iPts.push_back(1.);

    if (numIps > 1) {
        weights.push_back(2./(nIps*(nIps-1)) );
        for (int i=1; i<nIps-1;i++) {
            double lp = boost::math::legendre_p(nIps-1,iPts[i]);
            weights.push_back(2./(nIps*(nIps-1))/( lp * lp) );
        }
        weights.push_back(2./(nIps*(nIps-1)) );
    }
}

//! @brief returns the local coordinates of an integration point
//! @param rIpNum integration point (counting from zero)
//! @param rCoordinates (result)
Eigen::VectorXd NuTo::IntegrationType1D2NLobatto::GetLocalIntegrationPointCoordinates(int rIpNum) const
{
    if(rIpNum >= 0 && rIpNum < nIps)
        return Eigen::Matrix<double, 1, 1>::Constant(iPts[rIpNum]);
    else
        throw MechanicsException("[NuTo::IntegrationType1D2NLobatto::GetLocalIntegrationPointCoordinates] Ip number out of range.");
}


//! @brief returns the total number of integration points for this integration type
//! @return number of integration points
int NuTo::IntegrationType1D2NLobatto::GetNumIntegrationPoints()const
{
    return nIps;
}

//! @brief returns the weight of an integration point
//! @param rIpNum integration point (counting from zero)
//! @return weight of integration points
double NuTo::IntegrationType1D2NLobatto::GetIntegrationPointWeight(int rIpNum)const
{
    if(rIpNum >= 0 && rIpNum < nIps) return weights[rIpNum];
    throw MechanicsException("[NuTo::IntegrationType1D2NLobatto::GetIntegrationPointWeight] Ip number out of range.");
}

#ifdef ENABLE_VISUALIZE
void NuTo::IntegrationType1D2NLobatto::GetVisualizationCells(
    unsigned int& NumVisualizationPoints,
    std::vector<double>& VisualizationPointLocalCoordinates,
    unsigned int& NumVisualizationCells,
    std::vector<NuTo::eCellTypes>& VisualizationCellType,
    std::vector<unsigned int>& VisualizationCellsIncidence,
    std::vector<unsigned int>& VisualizationCellsIP) const
{
    NumVisualizationPoints = nIps+1;
    VisualizationPointLocalCoordinates.push_back(-1.);
    for (unsigned int i=1; i < NumVisualizationPoints-1; i++) {
        VisualizationPointLocalCoordinates.push_back(0.5 * (iPts[i-1] + iPts[i]));
    }
    VisualizationPointLocalCoordinates.push_back( 1.);

    NumVisualizationCells = nIps;
    for (unsigned int i=0; i < NumVisualizationCells; i++) {
        VisualizationCellType.push_back(NuTo::eCellTypes::LINE);
    }

    for (unsigned int i=0; i < NumVisualizationCells; i++) {
        VisualizationCellsIncidence.push_back(i);
        VisualizationCellsIncidence.push_back(i+1);
    }

    for (unsigned int i=0; i < NumVisualizationCells; i++) {
        VisualizationCellsIP.push_back(i);
    }
}
#endif // ENABLE_VISUALIZE

#ifdef ENABLE_SERIALIZATION
BOOST_CLASS_EXPORT_IMPLEMENT(NuTo::IntegrationType1D2NLobatto)
#endif
