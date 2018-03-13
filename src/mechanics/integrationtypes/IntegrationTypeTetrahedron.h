#pragma once

#include "mechanics/integrationtypes/IntegrationTypeBase.h"

namespace NuTo
{
//! @brief Integration types for the triangle
class IntegrationTypeTetrahedron : public IntegrationTypeBase
{

public:
    //! @brief data is taken from the attachment to the book "P. Solin, K. Segeth
    //! and I. Dolezel: Higher-Order Finite Element Methods", Chapman & Hall/CRC Press, 2003.
    //!
    //! They use another reference cell (running from -1 to 1) so the
    //! points and weights are changed accordingly
    static const std::vector<std::vector<Eigen::Vector4d>> quadratureData;

    static const std::vector<int> orderToIndex;

    //! @brief constructor
    //! @param order integration order
    IntegrationTypeTetrahedron(int order);

    //! @brief returns the local coordinates of an integration point
    //! @param i integration point (counting from zero)
    //! @return rCoordinates (result)
    Eigen::VectorXd GetLocalIntegrationPointCoordinates(int i) const override;

    //! @brief returns the total number of integration points for this integration type
    //! @return number of integration points
    int GetNumIntegrationPoints() const override;

    //! @brief returns the weight of an integration point
    //! @param i integration point (counting from zero)
    //! @return weight of integration points
    double GetIntegrationPointWeight(int i) const override;

    int GetDimension() const override
    {
        return 3;
    }

private:
    int mDataIndex;
};


} // namespace NuTo
