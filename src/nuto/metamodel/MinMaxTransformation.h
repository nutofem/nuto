// $Id$

/*******************************************************************************
Bauhaus-Universit�t Weimar
Author: Joerg F. Unger,  Septermber 2009
*******************************************************************************/


#ifndef MINMAXTRANSFORMATION_H
#define MINMAXTRANSFORMATION_H

#ifdef ENABLE_SERIALIZATION
#include <boost/serialization/access.hpp>
#include <boost/serialization/export.hpp>
#endif  // ENABLE_SERIALIZATION

#include "nuto/metamodel/Transformation.h"

namespace NuTo
{

//! @author Joerg F. Unger, ISM
//! @date September 2009
//! @brief abstract base class for the Transformations
class MinMaxTransformation : public Transformation
{
#ifdef ENABLE_SERIALIZATION
	friend class boost::serialization::access;
#endif  // ENABLE_SERIALIZATION

public:
	//! @brief constructor
	//! @param rCoordinate ... coordinate within the point coordinates (0<=entry<dim)
	//! @param rLb ... lower bound after the transformation
	//! @param rUb ... upper bound after the transformation
    MinMaxTransformation(unsigned int rCoordinate, double rLb, double rUb);
	
    //! @brief copy constructor
    //! @param rOther ... object which is copied
    MinMaxTransformation(const MinMaxTransformation &rOther);

    //! @brief destructor
	~MinMaxTransformation()
	{}

#ifdef ENABLE_SERIALIZATION
	//! @brief serializes the class
    //! @param ar         archive
    //! @param version    version
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version);
#endif  // ENABLE_SERIALIZATION

    //! @brief build the transformation using the given Points
    //! @param rCoordinates ... input point coordinates
    virtual void Build(const FullMatrix<double, Eigen::Dynamic, Eigen::Dynamic>& rCoordinates);

    //! @brief transform the given points in forward direction x = f(x) 
    //! @brief rCoordinates ... input point coordinates
    virtual void TransformForward(FullMatrix<double, Eigen::Dynamic, Eigen::Dynamic>& rCoordinates)const;

    //! @brief transform the given points in backward direction x = f^(-1)(x)
    //! @brief rCoordinates ... input point coordinates
    virtual void TransformBackward(FullMatrix<double, Eigen::Dynamic, Eigen::Dynamic>& rCoordinates)const;

protected:
    int  mCoordinate;     //!< coordinate within the point coordinates (0<=entry<dim
    double  mMin;        //!< min value of given coordinates
    double  mMax;        //!< max value of given coordinates
    double  mUb;         //!< upper bound after the transformation
    double  mLb;         //!< lower bound after the transformation

    // default constructor required by serialize
    MinMaxTransformation() : Transformation(){}
};
} // namespace nuto
#ifdef ENABLE_SERIALIZATION
#ifndef SWIG
BOOST_CLASS_EXPORT_KEY(NuTo::MinMaxTransformation)
#endif // SWIG
#endif  // ENABLE_SERIALIZATION

#endif /* MINMAXTRANSFORMATION_H */
