#ifdef ENABLE_SERIALIZATION
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/xml_oarchive.hpp>
#include <boost/archive/xml_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#endif // ENABLE_SERIALIZATION

#include "base/Exception.h"
#include "metamodel/MinMaxTransformation.h"


// constructor
NuTo::MinMaxTransformation::MinMaxTransformation(unsigned int rCoordinate, double rLb, double rUb) :
	NuTo::Transformation(),
	mCoordinate(rCoordinate),
	mMin(0.0),
	mMax(0.0),
	mUb(rUb),
	mLb(rLb)
{
}

NuTo::MinMaxTransformation::MinMaxTransformation(const MinMaxTransformation &rOther)
{
    mCoordinate = rOther.mCoordinate;
    mLb = rOther.mLb;
    mUb = rOther.mUb;
    mMin = rOther.mMin;
    mMax = rOther.mMax;
}


void NuTo::MinMaxTransformation::Build(const Eigen::MatrixXd& rCoordinates)
{
    if ( rCoordinates.cols() == 0)
	{
	    throw Exception("MinMaxTransformation::build - numberOfPoints must be greater than zero");
	}
    if ( rCoordinates.rows() <= mCoordinate)
    {
        throw Exception("MinMaxTransformation::build - coordinate to be transformed is out of range - check the dimension of your Matrix.");
    }
    const double *theptr = &rCoordinates.data()[mCoordinate];
	mMin = *theptr;
	mMax = *theptr;
    for (int count=1; count<rCoordinates.cols(); count++)
	{
        theptr+=rCoordinates.rows();
		if (*theptr<mMin)
			mMin = *theptr;
		else
		{
		    if (*theptr>mMax)
			    mMax = *theptr;
		}
	}
}

void NuTo::MinMaxTransformation::TransformForward(Eigen::MatrixXd& rCoordinates)const
{
    if ( rCoordinates.cols() == 0)
    {
        throw Exception("MinMaxTransformation::TransformForward - numberOfPoints must be greater than zero");
    }
    if ( rCoordinates.rows() <= mCoordinate)
    {
        throw Exception("MinMaxTransformation::TransformForward - coordinate to be transformed is out of range - check the dimension of your Matrix.");
    }
    double *theptr =  &rCoordinates.data()[mCoordinate];
    double deltaBound = mUb - mLb;
    double deltaValue = mMax-mMin;
	
	if (deltaBound==0 )
	{
        throw Exception("MinMaxTransformation::TransformForward - delta of prescribed bounds equal to zero");
	}

	if (deltaValue==0)
	{
        throw Exception("MinMaxTransformation::TransformForward - interval between min and max value of given points has size zero");
	}

    for (int count=0; count<rCoordinates.cols(); count++,theptr+=rCoordinates.rows())
	{
	    *theptr = mLb + (*theptr-mMin)/deltaValue*deltaBound;
	}
}

void NuTo::MinMaxTransformation::TransformBackward(Eigen::MatrixXd& rCoordinates)  const
{
    if ( rCoordinates.cols() == 0)
    {
        throw Exception("MinMaxTransformation::TransformBackward - numberOfPoints must be greater than zero");
    }
    if ( rCoordinates.rows() <= mCoordinate)
    {
        throw Exception("MinMaxTransformation::TransformBackward - coordinate to be transformed is out of range - check the dimension of your Matrix.");
    }
    double *theptr =  &rCoordinates.data()[mCoordinate];
    double deltaBound = mUb - mLb;
    double deltaValue = mMax-mMin;
	
	if (deltaBound==0 )
	{
        throw Exception("MinMaxTransformation::TransformBackward - delta of prescribed bounds equal to zero");
	}

	if (deltaValue==0)
	{
        throw Exception("MinMaxTransformation::TransformBackward - interval between min and max value of given points has size zero");
	}

    for (int count=0; count<rCoordinates.cols(); count++,theptr+=rCoordinates.rows())
	{
	    *theptr = mMin + (*theptr-mLb)/deltaBound*deltaValue;
	}
}

#ifdef ENABLE_SERIALIZATION
// serializes the class
template void NuTo::MinMaxTransformation::serialize(boost::archive::binary_oarchive & ar, const unsigned int version);
template void NuTo::MinMaxTransformation::serialize(boost::archive::xml_oarchive & ar, const unsigned int version);
template void NuTo::MinMaxTransformation::serialize(boost::archive::text_oarchive & ar, const unsigned int version);
template void NuTo::MinMaxTransformation::serialize(boost::archive::binary_iarchive & ar, const unsigned int version);
template void NuTo::MinMaxTransformation::serialize(boost::archive::xml_iarchive & ar, const unsigned int version);
template void NuTo::MinMaxTransformation::serialize(boost::archive::text_iarchive & ar, const unsigned int version);
template<class Archive>
void NuTo::MinMaxTransformation::serialize(Archive & ar, const unsigned int version)
{
#ifdef DEBUG_SERIALIZATION
    std::cout << "start serialize MinMaxTransformation" << std::endl;
#endif
    ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(Transformation)
       & BOOST_SERIALIZATION_NVP(mCoordinate)
       & BOOST_SERIALIZATION_NVP(mMin)
       & BOOST_SERIALIZATION_NVP(mMax)
       & BOOST_SERIALIZATION_NVP(mUb)
       & BOOST_SERIALIZATION_NVP(mLb);
#ifdef DEBUG_SERIALIZATION
    std::cout << "finish serialize MinMaxTransformation" << std::endl;
#endif
}
BOOST_CLASS_EXPORT_IMPLEMENT(NuTo::MinMaxTransformation)
#endif // ENABLE_SERIALIZATION
