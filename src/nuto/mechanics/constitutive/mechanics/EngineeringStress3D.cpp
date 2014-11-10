// $Id$
#include <iostream>

#ifdef ENABLE_SERIALIZATION
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/xml_oarchive.hpp>
#include <boost/archive/xml_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#endif // ENABLE_SERIALIZATION

#include <math.h>

#include "nuto/mechanics/constitutive/mechanics/EngineeringStress3D.h"

// constructor
NuTo::EngineeringStress3D::EngineeringStress3D(): ConstitutiveOutputBase::ConstitutiveOutputBase(), FullVector<double,6>()
{
    for (unsigned int count = 0; count < 6; count++)
    {
        (*this)[count] = 0.0;
    }
}

// number of components
unsigned int NuTo::EngineeringStress3D::GetNumberOfComponents() const
{
    return 6;
}

// get Engineering stress
const double* NuTo::EngineeringStress3D::GetData() const
{
    return data();
}

//! @brief ... sets the components of the Engineering stress tensor
//! @param ... components of the Engineering stress tensor (stored in an array)
//! @sa mEngineeringStress
void NuTo::EngineeringStress3D::SetData(double rData[6])
{
    for (unsigned int count = 0; count < 6; count++)
    {
    	(*this)[count] = rData[count];
    }
}


// info routine
void NuTo::EngineeringStress3D::Info(unsigned short rVerboseLevel) const
{
    std::cout << "    components of Engineering stress tensor (vector notation): "
              << (*this)[0] << ", " << (*this)[1] << ", " << (*this)[2] << ", "
              << (*this)[3] << ", " << (*this)[4] << ", " << (*this)[5] << std::endl;
}

double NuTo::EngineeringStress3D::Norm() const
{
	NuTo::FullVector<double,6> Stress;
	Stress(0) = (*this)[0];
	Stress(1) = (*this)[1];
	Stress(2) = (*this)[2];
	Stress(3) = (*this)[3];
	Stress(4) = (*this)[4];
	Stress(5) = (*this)[5];

    //*******************************************************************
    //*    NorM = second invariant = sqrt (0.5 * Stress_dev : Stress_dev)                                                *
    //*******************************************************************
    double invariante_2 = ((Stress(0)-Stress(1))*(Stress(0)-Stress(1))+
                            (Stress(0)-Stress(2))*(Stress(0)-Stress(2))+
                            (Stress(1)-Stress(2))*(Stress(1)-Stress(2)))/6.+
                             Stress(3)*Stress(3) + Stress(4)*Stress(4) + Stress(5)*Stress(5);
    return std::sqrt(invariante_2);
}

double NuTo::EngineeringStress3D::YieldSurfaceDruckerPrager3D(double rBeta, double rHP) const
{
	double invariante_1 = (*this)[0]+(*this)[1]+(*this)[2];

	double norm_2 = (*this).Norm();   /* norm_2 := J2 = sqrt(invariante_2) */

    //*******************************************************************
    //*    check yield first time                                            *
    //*******************************************************************
    double F_BETA = invariante_1 * rBeta/3. ;
    double F_FLOW = -rHP;

    //* calculate yield condition and yield condition flag (active or not)
    //* Drucker Prager
    return F_BETA + norm_2 + F_FLOW;
}

//! @brief calculates the first and second derivative of the Drucker Prager yield surface with respect to the stress
//! @param rdF_dSigma return value (first derivative)
//! @param rd2F_d2Sigma return value (second derivative)
//! @param rStress current stress
//! @param rBETA parameter of the Drucker Prager yield surface
//! @return false, if the stress is on the hydrostatic axis, otherwise true
//bool NuTo::ConstitutiveBase::YieldSurfaceDruckerPrager3DDerivatives(
//		NuTo::FullVector<double,6>& rdF_dSigma,
//		NuTo::FullMatrix<double,6,6>& rd2F_d2Sigma,
//		const NuTo::EngineeringStress3D& rStress, double rBETA) const
bool NuTo::EngineeringStress3D::YieldSurfaceDruckerPrager3DDerivatives(
		NuTo::FullVector<double,6>& rdF_dSigma,
		NuTo::FullMatrix<double,6,6>& rd2F_d2Sigma, double rBETA) const
{
	NuTo::FullVector<double,6> rdJ2_dSigma;   /* derivative of norm_2 := J2 with respect to rStress */

	double invariante_1 = (*this)[0]+(*this)[1]+(*this)[2];
	double norm_2 = (*this).Norm();
	double factor, func;
	const double eps = 1e-5;


	if (norm_2 > eps) {
		factor = 1./(6.*norm_2);
	} else {
		factor = (1./6.)*(1/eps)*(1.5 - 0.5*std::pow(norm_2/eps,2.));
		func = (1/eps)*(1.5 - 0.5*std::pow(norm_2/eps,2.));
	}

    /* Calculate derivative of norm_2:=J2 */
    rdJ2_dSigma(0) = factor * (2.*(*this)[0]-(*this)[1]-(*this)[2]);
    rdJ2_dSigma(1) = factor * (2.*(*this)[1]-(*this)[0]-(*this)[2]);
    rdJ2_dSigma(2) = factor * (2.*(*this)[2]-(*this)[0]-(*this)[1]);
    rdJ2_dSigma(3) = factor * (6.*(*this)[3]);
    rdJ2_dSigma(4) = factor * (6.*(*this)[4]);
    rdJ2_dSigma(5) = factor * (6.*(*this)[5]);

    /* Calculate derivative of the Drucker Prager yield surface */
    rdF_dSigma = rdJ2_dSigma;
    for (int i=0 ; i<3 ; i++)
    {
    	rdF_dSigma(i) += rBETA/3.;
    }

    /* Calculate second derivative of the Drucker Prager yield surface */
#define DELTA(i, j) ((i==j)?1:0)

    if (norm_2 > eps) {
    	factor = 1./(2.*norm_2*norm_2);
    	for ( int i = 0; i < 6; i++ )
    		for ( int j = 0; j < 6; j++ )
    		{
    			if (i < 3 && j < 3)
    			{
    				rd2F_d2Sigma(i,j) = norm_2*(DELTA(i, j) - 1./3.);
    				rd2F_d2Sigma(i,j) -= ((*this)[i] - invariante_1/3.)*rdJ2_dSigma(j);
    				rd2F_d2Sigma(i,j) *= factor;
    			}
    			if (i < 3 && j > 2)
    				{
    				rd2F_d2Sigma(i,j) = -factor*((*this)[i] - invariante_1/3.)*rdJ2_dSigma(j);
    				}
    			if (i > 2 && j < 3)
    				{
    				rd2F_d2Sigma(i,j) = -2.*factor*(*this)[i]*rdJ2_dSigma(j);
    				}
    			if (i > 2 && j > 2)
    				{
    				rd2F_d2Sigma(i,j) = norm_2*DELTA(i, j);
    				rd2F_d2Sigma(i,j) -= (*this)[i]*rdJ2_dSigma(j);
    				rd2F_d2Sigma(i,j) *= 2.*factor;
    				}
    		}	// for loop

	} else {
		factor = 1./(2.*eps*eps);
    	for ( int i = 0; i < 6; i++ )
    		for ( int j = 0; j < 6; j++ )
    		{
    			if (i < 3 && j < 3)
    			{
    				rd2F_d2Sigma(i,j) = eps*eps*func*(DELTA(i, j) - 1./3.);
    				rd2F_d2Sigma(i,j) -= ((*this)[i] - invariante_1/3.)*(norm_2/eps)*rdJ2_dSigma(j);
    				rd2F_d2Sigma(i,j) *= factor;
    			}
    			if (i < 3 && j > 2)
    				{
    				rd2F_d2Sigma(i,j) = -factor*((*this)[i] - invariante_1/3.)*(norm_2/eps)*rdJ2_dSigma(j);
    				}
    			if (i > 2 && j < 3)
    				{
    				rd2F_d2Sigma(i,j) = -2.*factor*(*this)[i]*(norm_2/eps)*rdJ2_dSigma(j);
    				}
    			if (i > 2 && j > 2)
    				{
    				rd2F_d2Sigma(i,j) = eps*eps*func*DELTA(i, j);
    				rd2F_d2Sigma(i,j) -= (*this)[i]*(norm_2/eps)*rdJ2_dSigma(j);
    				rd2F_d2Sigma(i,j) *= 2.*factor;
    				}
    		}	// for loop
	}	// end if

    return true;
}

#ifdef ENABLE_SERIALIZATION

//! @brief serializes the class
//! @param ar         archive
//! @param version    version
template void NuTo::EngineeringStress3D::serialize(boost::archive::binary_oarchive & ar, const unsigned int version);
template void NuTo::EngineeringStress3D::serialize(boost::archive::xml_oarchive & ar, const unsigned int version);
template void NuTo::EngineeringStress3D::serialize(boost::archive::text_oarchive & ar, const unsigned int version);
template void NuTo::EngineeringStress3D::serialize(boost::archive::binary_iarchive & ar, const unsigned int version);
template void NuTo::EngineeringStress3D::serialize(boost::archive::xml_iarchive & ar, const unsigned int version);
template void NuTo::EngineeringStress3D::serialize(boost::archive::text_iarchive & ar, const unsigned int version);
template<class Archive>
void NuTo::EngineeringStress3D::serialize(Archive & ar, const unsigned int version)
{
#ifdef DEBUG_SERIALIZATION
    std::cout << "start serialize EngineeringStress3D" << std::endl;
#endif
    ar &  BOOST_SERIALIZATION_BASE_OBJECT_NVP(ConstitutiveOutputBase);
    ar & boost::serialization::make_nvp ("EngineeringStrain3DEigen",boost::serialization::base_object< FullVector<double,6> > ( *this ));
#ifdef DEBUG_SERIALIZATION
    std::cout << "finish serialize EngineeringStress3D" << std::endl;
#endif
}
#endif // ENABLE_SERIALIZATION
