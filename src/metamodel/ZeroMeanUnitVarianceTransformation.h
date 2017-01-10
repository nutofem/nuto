#pragma once

#ifdef ENABLE_SERIALIZATION
#include <boost/serialization/access.hpp>
#include <boost/serialization/export.hpp>
#endif  // ENABLE_SERIALIZATION

// parent
#include "metamodel/Transformation.h"

namespace NuTo
{

//! @author Stefan Eckardt
//! @date February 2010
//! @brief zero mean, unit variance transformation
class ZeroMeanUnitVarianceTransformation : public Transformation
{
#ifdef ENABLE_SERIALIZATION
	friend class boost::serialization::access;
#endif  // ENABLE_SERIALIZATION

public:
    //! @brief constructor
    //! @param rCoordinate ... coordinate within the point coordinates
    //! @sa mCoordinate
    ZeroMeanUnitVarianceTransformation(unsigned int rCoordinate);
	
    //! @brief copy constructor
    //! @param other ... other object
    ZeroMeanUnitVarianceTransformation(const ZeroMeanUnitVarianceTransformation &other);

    //! @brief destructor
	~ZeroMeanUnitVarianceTransformation()
	{}

#ifdef ENABLE_SERIALIZATION
	//! @brief serializes the class
    //! @param ar         archive
    //! @param version    version
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version);
#endif  // ENABLE_SERIALIZATION

    //! @brief build the transformation using the given Points
    //! @param rCoordinates ... point coordinates
    virtual void Build(const Eigen::MatrixXd& rCoordinates);

    //! @brief transform the given points in forward direction x = f(x) 
    //! @param rCoordinates ... point coordinates
    virtual void TransformForward(Eigen::MatrixXd& rCoordinates)const;

    //! @brief transform the given points in backward direction x = f^(-1)(x)
    //! @param rCoordinates ... point coordinates
    virtual void TransformBackward(Eigen::MatrixXd& rCoordinates)const;

protected:
    int  mCoordinate;     //!< coordinate within the point coordinates (0<=entry<dim)
    double  mMean;        //!< mean value of given coordinates
    double  mStandardDeviation; //!< standard deviation of given coordinates

    //! @brief default constructor required by serialize
    ZeroMeanUnitVarianceTransformation(){}
};


} // namespace nuto
#ifdef ENABLE_SERIALIZATION
BOOST_CLASS_EXPORT_KEY(NuTo::ZeroMeanUnitVarianceTransformation)
#endif  // ENABLE_SERIALIZATION

