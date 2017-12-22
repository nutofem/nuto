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
    //! @param numIPs number of integration points
    //! @param method integration method (currently either GAUSS or LOBATTO)
    IntegrationTypeTensorProduct(size_t numIps, NuTo::eIntegrationMethod method);

    //! @brief returns the dimension
    int GetDimension() const override
    {
        return TDim;
    }

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

private:
    //! @brief ... 1D integration points coordinates
    std::vector<double> mIPts1D;
    //! @brief ... tensor product integration points coordinates
    std::vector<Eigen::Matrix<double, TDim, 1>> mIPts;
    //! @brief ... weights for the integration
    std::vector<double> mWeights;
};
} // namespace
