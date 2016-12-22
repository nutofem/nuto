// $Id$

#pragma once

#ifdef ENABLE_SERIALIZATION
#include <boost/serialization/access.hpp>
#include <boost/serialization/export.hpp>
#include <boost/serialization/vector.hpp>
#include <boost/serialization/array.hpp>
#else
#include <vector>
#endif //ENABLE_SERIALIZATION


namespace NuTo
{
#ifdef ENABLE_VISUALIZE
class CellBase;
enum class eCellTypes;
#endif // ENABLE_VISUALIZE



//! @author Daniel Arnold, ISM
//! @date Febrauary 2011
//! @brief ... standard abstract class for all modifyable integration points
class IntegrationPointBase
{
#ifdef ENABLE_SERIALIZATION
    friend class boost::serialization::access;
#endif  // ENABLE_SERIALIZATION

public:
#ifndef SWIG
#endif

    //! @brief constructor
    IntegrationPointBase();

    //! @brief constructor
    //! @param rCoords(Input):		Vector of local coordinates for this IP
    //! @param rWeight(Input):		This weighting factor
    //! @param rBoundingBox(Input):	Vector of the BoundingBox's coordinates [P1(1),P1(2),....,Pn(1),Pn(2),Pn(3)]
    IntegrationPointBase(const std::vector<double>& rCoords, const double& rWeight, const std::vector<double>& rBoundingBox);
    //! @brief destructor
    ~IntegrationPointBase();

    //! @brief Setter
    //! @param rWeight
    /**
     * Set the weight of this integration point
     */
    void SetWeight(const double &rWeight);

    //! @brief weight getter
    //! @return double weighting factor
    /**
     * return the weighting factor of this integration point
     */
    double GetWeight() const;
    //! @brief weight getter
    //! @return double weighting factor
    /**
     * return the weighting factor of this integration point
     */
    void GetWeight(double& rWeight);

    //! @brief coordinate getter
    //! @return vector<double> coordinate vector
    /**
     * return the coordinates of this integration point
     */
    std::vector<double> GetLocalCoords() const;

    //! @brief coordinate getter
    //! @param vector<double> coordinate vector
    /**
     * return the coordinates of this integration point
     */
    void GetLocalCoords(std::vector<double>& rCoords);

#ifdef ENABLE_VISUALIZE
    //! @brief getter of the visualization cell
    //! @param double weighting factor
    /**
     * return the information of this integration cell:
     * contains the information of the number of belonging visualization points, the cell type,
     * the coordiantes and the incidences
     */
    void GetVisualizationCell(unsigned int rNumVisualizationPoints,
                                eCellTypes rVisualizationCellType,
                                std::vector<double>& rVisualizationPointLocalCoordinates,
                                std::vector<unsigned int>& rVisualizationCellsIncidence ) const;
#endif // ENABLE_VISUALIZE


#ifdef ENABLE_SERIALIZATION
    //! @brief serializes the class
    //! @param ar         archive
    //! @param version    version
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version);
#endif // ENABLE_SERIALIZATION


protected:

    double mWeight;
    std::vector<double> mCoords;
#ifdef ENABLE_VISUALIZE
    unsigned int mNumVisualizationPoints;
    NuTo::eCellTypes mVisualizationCellType;
    std::vector<double> mVisualizationPointLocalCoordinates;
    std::vector<unsigned int> mVisualizationCellsIncidence;
#endif // ENABLE_VISUALIZE

};

}//namespace NuTo
#ifdef ENABLE_SERIALIZATION
BOOST_CLASS_EXPORT_KEY(NuTo::IntegrationPointBase)
BOOST_CLASS_TRACKING(NuTo::IntegrationPointBase, track_always)
#endif // ENABLE_SERIALIZATION

