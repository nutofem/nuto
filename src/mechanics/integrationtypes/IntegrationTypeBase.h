#pragma once

#ifdef ENABLE_SERIALIZATION
#include <boost/serialization/access.hpp>
#include <boost/serialization/export.hpp>
#include <iostream>
#endif // ENABLE_SERIALIZATION

#include <vector>
#include <eigen3/Eigen/Core>
#include "mechanics/MechanicsException.h"

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
#ifdef ENABLE_SERIALIZATION
    friend class boost::serialization::access;
#endif // ENABLE_SERIALIZATION

public:
    IntegrationTypeBase() = default;
    IntegrationTypeBase(const IntegrationTypeBase&) = default;
    IntegrationTypeBase(IntegrationTypeBase&&) = default;

    IntegrationTypeBase& operator=(const IntegrationTypeBase&) = default;
    IntegrationTypeBase& operator=(IntegrationTypeBase&&) = default;

    //! @brief ... destructor
    virtual ~IntegrationTypeBase() = default;

#ifdef ENABLE_SERIALIZATION
    //! @brief serializes the class
    //! @param ar         archive
    //! @param version    version
    template <class Archive>
    void serialize(Archive& ar, const unsigned int version)
    {
#ifdef DEBUG_SERIALIZATION
        std::cout << "start serialize IntegrationTypeBase" << std::endl;
#endif
#ifdef DEBUG_SERIALIZATION
        std::cout << "finish serialize IntegrationTypeBase" << std::endl;
#endif
    }
#endif // ENABLE_SERIALIZATION

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

    virtual void GetVisualizationCells(unsigned int& NumVisualizationPoints,
                                       std::vector<double>& VisualizationPointLocalCoordinates,
                                       unsigned int& NumVisualizationCells,
                                       std::vector<NuTo::eCellTypes>& VisualizationCellType,
                                       std::vector<unsigned int>& VisualizationCellsIncidence,
                                       std::vector<unsigned int>& VisualizationCellsIP) const {};
#endif // ENABLE_VISUALIZE
protected:
};
} // namespace NuTo

#ifdef ENABLE_SERIALIZATION
BOOST_CLASS_EXPORT_KEY(NuTo::IntegrationTypeBase)
#endif // ENABLE_SERIALIZATION
