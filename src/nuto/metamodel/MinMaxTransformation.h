/*******************************************************************************
Bauhaus-Universit�t Weimar
Author: Joerg F. Unger,  Septermber 2009
*******************************************************************************/


#ifndef MINMAXTRANSFORMATION_H
#define MINMAXTRANSFORMATION_H

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

    MinMaxTransformation(unsigned int rCoordinate, double rLb, double rUb) : Transformation()
	{
	    mCoordinate = rCoordinate;
	    mLb = rLb;
	    mUb = rUb;
	}
	
    MinMaxTransformation(const MinMaxTransformation &other)
	{
	    mCoordinate = other.mCoordinate;
	    mLb = other.mLb;
	    mUb = other.mUb;
	}
    
	~MinMaxTransformation()
	{}

#ifdef ENABLE_SERIALIZATION
	//! @brief serializes the class
    //! @param ar         archive
    //! @param version    version
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version)
    {
        ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(Transformation)
           & BOOST_SERIALIZATION_NVP(mCoordinate)
           & BOOST_SERIALIZATION_NVP(mMin)
           & BOOST_SERIALIZATION_NVP(mMax)
           & BOOST_SERIALIZATION_NVP(mUb)
           & BOOST_SERIALIZATION_NVP(mLb);
    }
#endif  // ENABLE_SERIALIZATION

    //! @brief build the transformation using the given Points
    virtual void Build(const FullMatrix<double>& rCoordinates);

    //! @brief transform the given points in forward direction x = f(x) 
    virtual void TransformForward(FullMatrix<double>& rCoordinates)const;

    //! @brief transform the given points in backward direction x = f^(-1)(x)
    virtual void TransformBackward(FullMatrix<double>& rCoordinates)const;

protected:
    int  mCoordinate;     //!< coordinate within the point coordinates (0<=entry<dim
    double  mMin;        //!< min value of given coordinates
    double  mMax;        //!< max value of given coordinates
    double  mUb;         //!< upper bound after the transformation
    double  mLb;         //!< lower bound after the transformation
};


} // namespace nuto

#endif /* MINMAXTRANSFORMATION_H */
