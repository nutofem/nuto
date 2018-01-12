#pragma once

#include <vector>
#include "nuto/mechanics/integrationtypes/IntegrationTypeBase.h"
#include "nuto/mechanics/integrationtypes/IntegrationTypeEnum.h"
#include "nuto/mechanics/elements/ElementEnum.h"

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
    int GetCoordinateDimension() const override
    {
        return TDim;
    }

    //! @brief returns the local coordinates of an integration point
    //! @param rIpNum integration point (counting from zero)
    //! @return rCoordinates (result)
    Eigen::VectorXd GetLocalIntegrationPointCoordinates(int rIpNum) const;

    //! @brief returns the total number of integration points for this integration type
    //! @return number of integration points
    int GetNumIntegrationPoints() const override;

    std::string GetStrIdentifier() const override
    {
        std::string identifier;
        if (TDim == 2)
        {
            if (mMethod == NuTo::eIntegrationMethod::GAUSS)
                identifier =
                        std::string("2D4NGAUSS") + std::to_string(mIPts1D.size() * mIPts1D.size()) + std::string("IP");
            else
                identifier = std::string("2D4NLOBATTO") + std::to_string(mIPts1D.size() * mIPts1D.size()) +
                             std::string("IP");
        }

        return identifier;
    }

    virtual bool CheckElementCompatibility(NuTo::Element::eElementType rElementType) const override
    {
        switch (rElementType)
        {
        case NuTo::Element::eElementType::CONTINUUMELEMENT:
        case NuTo::Element::eElementType::CONTINUUMBOUNDARYELEMENT:
        case NuTo::Element::eElementType::CONTINUUMCONTACTELEMENT:
        case NuTo::Element::eElementType::CONTINUUMBOUNDARYELEMENTCONSTRAINEDCONTROLNODE:
            return true;
        default:
            return false;
        }
    }

    //! @brief returns the local coordinates of an integration point
    //! @param rIpNum integration point (counting from zero)
    //! @param rCoordinates (result)
    void GetLocalIntegrationPointCoordinates2D(int rIpNum, double rCoordinates[2]) const
    {
        if (rIpNum >= 0 && rIpNum < (int)mIPts.size())
        {
            rCoordinates[0] = mIPts[rIpNum](0, 0);
            rCoordinates[1] = mIPts[rIpNum](1, 0);
        }
        else
            throw Exception(__PRETTY_FUNCTION__, "Ip number out of range.");
    }

    //! @brief returns the weight of an integration point
    //! @param rIpNum integration point (counting from zero)
    //! @return weight of integration points
    double GetIntegrationPointWeight(int rIpNum) const override;

    void GetVisualizationCells(unsigned int& NumVisualizationPoints,
                               std::vector<double>& VisualizationPointLocalCoordinates,
                               unsigned int& NumVisualizationCells,
                               std::vector<NuTo::eCellTypes>& VisualizationCellType,
                               std::vector<unsigned int>& VisualizationCellsIncidence,
                               std::vector<unsigned int>& VisualizationCellsIP) const override;

private:
    //! @brief ... 1D integration points coordinates
    std::vector<double> mIPts1D;
    //! @brief ... tensor product integration points coordinates
    std::vector<Eigen::Matrix<double, TDim, 1>> mIPts;
    //! @brief ... weights for the integration
    std::vector<double> mWeights;
    //! @brief ... method
    NuTo::eIntegrationMethod mMethod;

    //! @brief ... visualization points in 1D
    void GetVisualizationPoints(unsigned int& NumVisualizationPoints,
                                std::vector<double>& VisualizationPointLocalCoordinates) const;
};
} // namespace
