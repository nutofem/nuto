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
#include "math/Legendre.h"

// constructor
NuTo::IntegrationType1D2NLobatto::IntegrationType1D2NLobatto(int nIps)
{
    mIPts = NuTo::Math::Polynomial::LegendreDerivRoots(nIps-1);
    mIPts.insert(mIPts.begin(),-1.);
    mIPts.push_back(1.);

    if (nIps > 1) {
        mWeights.push_back(2./(nIps*(nIps-1)) );
        for (int i=1; i<nIps-1;i++) {
            double lp = NuTo::Math::Polynomial::Legendre(nIps-1,mIPts[i]);
            mWeights.push_back(2./(nIps*(nIps-1))/( lp * lp) );
        }
        mWeights.push_back(2./(nIps*(nIps-1)) );
    }
}

Eigen::VectorXd NuTo::IntegrationType1D2NLobatto::GetLocalIntegrationPointCoordinates(int rIpNum) const
{
    if(rIpNum >= 0 && rIpNum < mIPts.size())
        return Eigen::Matrix<double, 1, 1>::Constant(mIPts[rIpNum]);
    else
        throw MechanicsException("[NuTo::IntegrationType1D2NLobatto::GetLocalIntegrationPointCoordinates] Ip number out of range.");
}

int NuTo::IntegrationType1D2NLobatto::GetNumIntegrationPoints()const
{
    return mIPts.size();
}

double NuTo::IntegrationType1D2NLobatto::GetIntegrationPointWeight(int rIpNum)const
{
    if(rIpNum >= 0 && rIpNum < mIPts.size()) return mWeights[rIpNum];
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
    NumVisualizationPoints = mIPts.size()+1;
    VisualizationPointLocalCoordinates.push_back(-1.);
    for (unsigned int i=1; i < NumVisualizationPoints-1; i++) {
        VisualizationPointLocalCoordinates.push_back(0.5 * (mIPts[i-1] + mIPts[i]));
    }
    VisualizationPointLocalCoordinates.push_back( 1.);

    NumVisualizationCells = mIPts.size();
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
