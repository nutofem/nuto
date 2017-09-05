#pragma once

#include <vector>
#include <eigen3/Eigen/Core>
#include "base/Exception.h"

#ifdef ENABLE_VISUALIZE
#include "visualize/VisualizeEnum.h"
#endif // ENABLE_VISUALIZE


namespace NuTo
{

class IntegrationPointBase;

namespace Element
{
enum class eElementType;
} // namespace Element


//! @author Jörg F. Unger, ISM
//! @date November 2009
//! @brief ... standard abstract class for all integration types
class IntegrationTypeBase
{

public:
    IntegrationTypeBase() = default;
    IntegrationTypeBase(const IntegrationTypeBase&) = default;
    IntegrationTypeBase(IntegrationTypeBase&&) = default;

    IntegrationTypeBase& operator=(const IntegrationTypeBase&) = default;
    IntegrationTypeBase& operator=(IntegrationTypeBase&&) = default;

    //! @brief ... destructor
    virtual ~IntegrationTypeBase() = default;

    virtual int GetDimension() const = 0;

    //! @brief returns the local coordinates of an integration point
    //! @param rIpNum integration point (counting from zero)
    //! @return rCoordinates (result)
    virtual Eigen::VectorXd GetLocalIntegrationPointCoordinates(int rIpNum) const = 0;


    virtual Eigen::MatrixXd GetNaturalIntegrationPointCoordinates() const;

    //! @brief returns the total number of integration points for this integration type
    //! @return number of integration points
    virtual int GetNumIntegrationPoints() const = 0;

    //! @brief returns the weight of an integration point
    //! @param rIpNum integration point (counting from zero)
    //! @return weight of integration points
    virtual double GetIntegrationPointWeight(int rIpNum) const = 0;

    //! @brief info about the integration type
    //! @param rVerboseLevel determines how detailed the information is
    void Info(int rVerboseLevel) const;


#ifdef ENABLE_VISUALIZE

    struct CellInfo
    {
        int visualizeCellId = -1; // global visualize id, set later by visualizer
        std::vector<int> pointIds;
        eCellTypes cellType;
        int ipId;
    };

    struct CellVectexInfo
    {
        int visualizePointId = -1; // global visualize id, set later by visualizer
        Eigen::VectorXd localCoords;
    };

    struct IpCellInfo
    {
        std::vector<CellInfo> cells;
        std::vector<CellVectexInfo> vertices;
    };

    virtual IpCellInfo GetVisualizationCells() const;

    virtual void GetVisualizationCells(unsigned int&, std::vector<double>&, unsigned int&,
                                       std::vector<NuTo::eCellTypes>&, std::vector<unsigned int>&,
                                       std::vector<unsigned int>&) const {};
#endif // ENABLE_VISUALIZE
protected:
};
} // namespace NuTo
