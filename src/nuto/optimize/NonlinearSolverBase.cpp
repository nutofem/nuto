// $Id$

#ifdef ENABLE_SERIALIZATION
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/xml_oarchive.hpp>
#include <boost/archive/xml_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <boost/serialization/vector.hpp>
#include <boost/ptr_container/serialize_ptr_map.hpp>
#include <boost/ptr_container/serialize_ptr_list.hpp>
#endif // ENABLE_SERIALIZATION


#include "nuto/optimize/NonlinearSolverBase.h"
#include "nuto/optimize/OptimizeException.h"

#include "nuto/math/FullVector.h"
#include "nuto/math/FullMatrix.h"

using namespace std;

NuTo::NonlinearSolverBase::NonlinearSolverBase() : NuTo::NuToObject::NuToObject()
{
	mTolResidual = 1.0e-12;
	mTolSolution = numeric_limits<double>::epsilon();
	mMaxIterationsNumber = 100;
	mResidualFunction = 0;
	mAssignResidual = false;
}

//! @brief ... numerical differentiation of the residual function mResidualFunction (numerical Jacobi matrix)
//! @param rParam ... parameters necessary to evaluate the residual
//! @param rUnknown ... position at which the derivative is taken
//! @param rFvec ... residual vector at the position rUnknown, rFvec = mResidualFunction(rParameter,rUnknown)
NuTo::FullMatrix<double,Eigen::Dynamic,Eigen::Dynamic> NuTo::NonlinearSolverBase::DResidualNum(NuTo::FullVector<double,Eigen::Dynamic> rUnknown,
		NuTo::FullVector<double,Eigen::Dynamic> &rFvec) const
{
	const double EPS = 1.0e-8;
	int n=rUnknown.GetNumRows();
	NuTo::FullMatrix<double,Eigen::Dynamic,Eigen::Dynamic> deriv(n,n);
	NuTo::FullVector<double,Eigen::Dynamic> xh=rUnknown;

	if (mResidualFunction==0 && mAssignResidual==false) {
		throw OptimizeException("[NuTo::NonLinearSolverBase::DResidualNum] the pointer to the residual function is required.");
	}

	for (int j=0;j<n;j++) {
		double temp=xh[j];
		double h=EPS*std::abs(temp);
	   	if (h == 0.0) h=EPS;
	   		xh[j]=temp+h;
	   		h=xh[j]-temp;
//  	   		NuTo::FullVector<double,Eigen::Dynamic> f=(*mResidualFunction)(this->mParameter,xh);
  	   		NuTo::FullVector<double,Eigen::Dynamic> f=(mResidualFunctionBoost)(this->mParameter,xh);
  	   		xh[j]=temp;
  	   		for (int i=0;i<n;i++)
 	   			deriv(i,j)=(f[i]-rFvec[i])/h;
	}
   	return deriv;
}

//! @brief ... calculates 0.5*rFvec^2 and updates rFvec = mResidualFunction(rParameter,rUnknown)
double NuTo::NonlinearSolverBase::Fmin(NuTo::FullVector<double,Eigen::Dynamic> rUnknown, NuTo::FullVector<double,Eigen::Dynamic> &rFvec) const
{
	if (mResidualFunction==0 && mAssignResidual==false) {
		throw OptimizeException("[NuTo::NonLinearSolverBase::DResidualNum] the pointer to the residual function is required.");
	}

//	rFvec = (*mResidualFunction)(this->mParameter,rUnknown);
	rFvec = (mResidualFunctionBoost)(this->mParameter,rUnknown);
	double sum=0;
	sum = rFvec.dot(rFvec);
	return 0.5*sum;
}

//! @brief ... Info routine that prints general information about the object (detail according to verbose level)
void NuTo::NonlinearSolverBase::Info()const
{
}

#ifdef ENABLE_SERIALIZATION
template void NuTo::NonlinearSolverBase::serialize(boost::archive::binary_oarchive & ar, const unsigned int version);
template void NuTo::NonlinearSolverBase::serialize(boost::archive::xml_oarchive & ar, const unsigned int version);
template void NuTo::NonlinearSolverBase::serialize(boost::archive::text_oarchive & ar, const unsigned int version);
template void NuTo::NonlinearSolverBase::serialize(boost::archive::binary_iarchive & ar, const unsigned int version);
template void NuTo::NonlinearSolverBase::serialize(boost::archive::xml_iarchive & ar, const unsigned int version);
template void NuTo::NonlinearSolverBase::serialize(boost::archive::text_iarchive & ar, const unsigned int version);
template<class Archive>
void NuTo::NonlinearSolverBase::serialize(Archive & ar, const unsigned int version)
{
#ifdef DEBUG_SERIALIZATION
    std::cout << "start serialization of NonlinearSolverBase" << "\n";
#endif
    ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(NuToObject)
       & BOOST_SERIALIZATION_NVP(mTolResidual)
       & BOOST_SERIALIZATION_NVP(mTolSolution)
       & BOOST_SERIALIZATION_NVP(mMaxIterationsNumber)
       & BOOST_SERIALIZATION_NVP(mParameter);
//       & BOOST_SERIALIZATION_NVP(mResidualFunction);
#ifdef DEBUG_SERIALIZATION
    std::cout << "finish serialization of NonlinearSolverBase" << "\n";
#endif
}
#endif  // ENABLE_SERIALIZATION

