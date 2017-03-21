#ifdef ENABLE_SERIALIZATION
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/xml_oarchive.hpp>
#include <boost/archive/xml_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <boost/serialization/vector.hpp>
#include <boost/ptr_container/serialize_ptr_map.hpp>
#endif // ENABLE_SERIALIZATION

# ifdef _OPENMP
#include <omp.h>
# endif

#include "optimize/NewtonRaphson.h"
#include "optimize/OptimizeException.h"

#include <eigen3/Eigen/LU>

//! @brief constructor
NuTo::NewtonRaphson::NewtonRaphson() : NonlinearSolverBase()
{
	mCheckNewtonRaphson = false;
    mResidualDerivativeFunction = nullptr;
    mAssignResidualResidualDerivative = false;
}


//#ifdef ENABLE_SERIALIZATION
//// serializes the class
//template void NuTo::NewtonRaphson::serialize(boost::archive::binary_oarchive & ar, const unsigned int version);
//template void NuTo::NewtonRaphson::serialize(boost::archive::xml_oarchive & ar, const unsigned int version);
//template void NuTo::NewtonRaphson::serialize(boost::archive::text_oarchive & ar, const unsigned int version);
//template void NuTo::NewtonRaphson::serialize(boost::archive::binary_iarchive & ar, const unsigned int version);
//template void NuTo::NewtonRaphson::serialize(boost::archive::xml_iarchive & ar, const unsigned int version);
//template void NuTo::NewtonRaphson::serialize(boost::archive::text_iarchive & ar, const unsigned int version);
//template<class Archive>
//void NuTo::NewtonRaphson::serialize(Archive & ar, const unsigned int version)
//{
//    #ifdef DEBUG_SERIALIZATION
//        std::cout << "start serialization of NewtonRaphson" << "\n";
//    #endif
//        ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(NonlinearSolverBase);
//           & BOOST_SERIALIZATION_NVP(mResidualDerivativeFunction);
//    #ifdef DEBUG_SERIALIZATION
//        std::cout << "finish serialization of NewtonRaphson" << "\n";
//    #endif
//}
//#endif // ENABLE_SERIALIZATION


void NuTo::NewtonRaphson::Solve(Eigen::VectorXd &rUnknown)
{
	if (not mAssignResidual && mResidualFunction == nullptr) {
		throw OptimizeException("[NuTo::NewtonRaphson::Solve] the pointer to the residual function is required.");
	}

	this->NewtonRaphsonIterator(rUnknown,this->mCheckNewtonRaphson);
}

//! @brief ... the routine performs Newton-Raphson iterations
void NuTo::NewtonRaphson::NewtonRaphsonIterator(Eigen::VectorXd &rX, bool &rCheck) const
{
	const int MAXITS(this->mMaxIterationsNumber);
	const double TOLF = this->mTolResidual, TOLMIN=1.0e-12, STPMX=100.0;
	const double TOLX = this->mTolSolution;
	int its,n=rX.rows();
	double den,f,fold,stpmax,test;
	Eigen::VectorXd g(n),p(n),xold(n);
    Eigen::MatrixXd fjac(n,n);
	Eigen::VectorXd fvec;

	f = this->Fmin(rX, fvec);

	test = fvec.array().abs().maxCoeff();

//	std::cout << "	fvec_norm before LineSearch " << test << std::endl;	// Tests

	if (test < 0.01*TOLF)
    {
		rCheck=false;
        return;
	}

	stpmax=STPMX*std::max(rX.lpNorm<2>(),double(n));

	for (its=0;its<MAXITS;its++) {
//   		std::cout<< "===== Iteration =====" << its <<std::endl;   // Test
		if (this->mResidualDerivativeFunction != nullptr || this->mAssignResidualResidualDerivative == true) {
			// if analytical Jacobi is given
//			fjac = (*mResidualDerivativeFunction)(this->mParameter,rX);
			fjac = (mResidualDerivativeFunctionBoost)(this->mParameter,rX);
//   			std::cout<<"*** Analytical ***"<<std::endl;  	// Test
//   			std::cout << fjac << std::endl;             	// Test
		} else {
			// takes numerical Jacobi, if analytical Jacobi is not given
			fjac=this->DResidualNum(rX,fvec);
		}

		g = fjac.transpose() * fvec;

		xold = rX;
		fold=f;
		p = -fvec;

// 		std::cout << "	fvec before LineSearch = " << fvec.transpose()<< std::endl;   // Test
// 		std::cout << "	p before LineSearch    = " << fjac.fullPivLu().solve(-fvec).transpose()<< std::endl;   // Test

		p =	fjac.fullPivLu().solve(-fvec).transpose();		// LU SOLVER of fjac * p = -fvec
//															// SVD SOLVER
//		p = fjac.jacobiSvd().solve(-fvec).transpose();      // SVD SOLVER

		this->LineSearch(xold,fold,g,p,rX,f,stpmax,rCheck,fvec);

//      std::cout << "	fvec_norm after LineSearch " << fvec.array().abs().maxCoeff() << std::endl;	 // Test
//    	std::cout << "	fvec after LineSearch = " << fvec.transpose()<<std::endl;  					 // Test
//    	std::cout << "	x    after LineSearch = " << rX.transpose()<<std::endl;  					 // Test

		test = fvec.array().abs().maxCoeff();

		if (test < 0.01*TOLF)
        {
			// ordinary return
			rCheck=false;
            return;
		}

		if (rCheck)
        {
			// spurious solution, local minimum of Fmin
			den=std::max(f,0.5*n);
			test = ( g.array().abs() * rX.array().abs().max(1.0) ).maxCoeff() / den;
			rCheck=(test < TOLMIN);
            return;
		}

		test = ( (rX.array()-xold.array()).abs() / rX.array().abs().max(1.0) ).maxCoeff();   // OPTIMIZED

		if (test < TOLX)
        {
			// too small change of the solution, spurious solution
            return;
        }
	}
    throw NuTo::OptimizeException(__PRETTY_FUNCTION__, "The maximal number of iterations exceeded");
}

//! @brief ... the routine performs line search correction of the Newton step
void NuTo::NewtonRaphson::LineSearch(Eigen::VectorXd &rXold, const double rFold, Eigen::VectorXd &rG,
		Eigen::VectorXd &rP, Eigen::VectorXd &rX, double &rF, const double rStpmax, bool &rCheck,
		Eigen::VectorXd &rFvec) const
{
	const double ALF=1.0e-4, TOLX = this->mTolSolution;
	double a,alam,alam2=0.0,alamin,b,disc,f2=0.0;
	double rhs1,rhs2,slope=0.0,sum=0.0,test,tmplam;

	rCheck = false;

	sum = rP.lpNorm<2>();

	if (sum > rStpmax)
		rP *= rStpmax/sum;

	slope = rG.dot(rP);

	if (slope >= 0.0){
		std::cout << "Roundoff problem in linesearch." << std::endl;
		throw("[NuTo::NewmarkRaphson::LineSearch] Roundoff problem in linesearch.");
	}

	test=0.0;

	test = ( rP.array().abs() / rXold.array().abs().max(1.0) ).maxCoeff();

	alamin=TOLX/test;
	alam=1.0;
	for (;;) {
		rX = rXold + alam*rP;
		rF = this->Fmin(rX,rFvec);

		if (alam < alamin) {

			rX = rXold;
			rCheck=true;
			return;
		} else if (rF <= rFold+ALF*alam*slope) return;
		else {
			if (alam == 1.0)
				tmplam = -slope/(2.0*(rF-rFold-slope));
			else {
				rhs1=rF-rFold-alam*slope;
				rhs2=f2-rFold-alam2*slope;
				a=(rhs1/(alam*alam)-rhs2/(alam2*alam2))/(alam-alam2);
				b=(-alam2*rhs1/(alam*alam)+alam*rhs2/(alam2*alam2))/(alam-alam2);
				if (a == 0.0) tmplam = -slope/(2.0*b);
				else {
					disc=b*b-3.0*a*slope;
					if (disc < 0.0) tmplam=0.5*alam;
					else if (b <= 0.0) tmplam=(-b+sqrt(disc))/(3.0*a);
					else tmplam=-slope/(b+sqrt(disc));
				}
				if (tmplam>0.5*alam)
					tmplam=0.5*alam;
			}
		}
		alam2=alam;
		f2 = rF;
		alam=std::max(tmplam,0.1*alam);
	}
}

#ifdef ENABLE_SERIALIZATION
#ifndef SWIG
BOOST_CLASS_EXPORT_IMPLEMENT(NuTo::NewtonRaphson)
#endif // SWIG
#endif // ENABLE_SERIALIZATION
