#pragma once

#include <vector>
#include "mechanics/integrationtypes/IntegrationTypeBase.h"
#include "mechanics/integrationtypes/IntegrationTypeEnum.h"

namespace NuTo
{
//! @author Philipp Müller, Peter Otto, BAM
//! @date Jul 2017
//! @brief ... integration types in 1,2,3D; tensor product of 1D Lobatto or Gauss

template <int TDim>
class IntegrationTypeTensorProduct : public IntegrationTypeBase
{
public:
    //! @brief constructor
    IntegrationTypeTensorProduct(size_t numIps, NuTo::eIntegrationMethod method);

    //! @brief computes points and weights for Lobatto quadrature in 1D
    //! @param numIPs number of integration points
    //! @return pair of quadrature weights and points range [-1,1] including boundary points
    std::pair<std::vector<double>, std::vector<double>>  ComputeWeightsAndPoints1DLobatto(int nIps);

    //! @brief computes points and weights for Gauss quadrature in 1D
    //! @param numIPs number of integration points
    //! @return pair of quadrature weights and points range (-1,1)
    std::pair<std::vector<double>, std::vector<double> > ComputeWeightsAndPoints1DGauss(int nIps);

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
    //! @brief ... 1D integration points coordinates
    std::vector<double> mIPts1D;
    //! @brief ... tensor product integration points coordinates
    std::vector<Eigen::Matrix<double, TDim, 1>> mIPts;
    //! @brief ... weights for the integration
    std::vector<double> mWeights;

#ifdef ENABLE_VISUALIZE
    //! @brief ... visualization points in 1D
    void GetVisualizationPoints(unsigned int& NumVisualizationPoints, std::vector<double>& VisualizationPointLocalCoordinates) const;
#endif // ENABLE_VISUALIZE
};
} // namespace
