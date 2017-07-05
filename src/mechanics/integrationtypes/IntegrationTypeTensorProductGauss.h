#pragma once

#include <vector>
#include "mechanics/integrationtypes/IntegrationType1D2NGauss.h"

namespace NuTo
{
//! @author Philipp Müller, Peter Otto, BAM
//! @date Jul 2017
//! @brief ... integration types in 1,2,3D with Gauss integration, tensor product of 1D
template <int TDim>
class IntegrationTypeTensorProductGauss : public IntegrationTypeBase
{
public:
    //! @brief constructor
    IntegrationTypeTensorProductGauss(int numIps);


    //! @brief returns the dimension
    int GetDimension() const override {return TDim;}

    //! @brief returns the local coordinates of an integration point
    //! @param rIpNum integration point (counting from zero)
    //! @return rCoordinates (result)
    Eigen::VectorXd GetLocalIntegrationPointCoordinates(int rIpNum) const override;

    //! @brief returns the total number of integration points for this integration type
    //! @return number of integration points
    int GetNumIntegrationPoints() const override;

    //! @brief returns the weight of an integration point
    //! @param rIpNum integration point (counting from zero)
    //! @return weight of integration points
    double GetIntegrationPointWeight(int rIpNum) const override;

#ifdef ENABLE_VISUALIZE
    void GetVisualizationCells(unsigned int& NumVisualizationPoints,
                               std::vector<double>& VisualizationPointLocalCoordinates,
                               unsigned int& NumVisualizationCells,
                               std::vector<NuTo::eCellTypes>& VisualizationCellType,
                               std::vector<unsigned int>& VisualizationCellsIncidence,
                               std::vector<unsigned int>& VisualizationCellsIP) const override;
#endif // ENABLE_VISUALIZE
private:
    //! @brief ... integration points coordinates
    std::vector<Eigen::Matrix<double, TDim, 1>> mIPts;
    //! @brief ... weights for the integration
    std::vector<double> mWeights;
    //! @brief ... number of IPs in 1D
    int mNumIPs1D;
};
} // namespace
