// $Id$

#ifndef DAMAGEVISCOPLASTICITYENGINEERINGSTRESS_H_
#define DAMAGEVISCOPLASTICITYENGINEERINGSTRESS_H_

#include "nuto/mechanics/constitutive/ConstitutiveEnum.h"
#include "nuto/mechanics/elements/ElementEnum.h"
#include "nuto/mechanics/constitutive/ConstitutiveBase.h"
#include "nuto/mechanics/constitutive/ConstitutiveTangentLocal.h"
//#include "nuto/optimize/NewtonRaphson.h"

using namespace std;

namespace NuTo
{
//! @brief ... visco-plastic material model with damage
class ConstitutiveStaticDataDamageViscoPlasticity3D;
class Logger;
template <int TNumRows,int TNumColumns>
class ConstitutiveTangentLocal;

class DamageViscoPlasticityEngineeringStress: public ConstitutiveBase
{
#ifdef ENABLE_SERIALIZATION
    friend class boost::serialization::access;
#endif // ENABLE_SERIALIZATION
public:
    DamageViscoPlasticityEngineeringStress();

    //! @brief ... evaluate the constitutive relation in 2D
    //! @param rElement ... element
    //! @param rIp ... integration point
    //! @param rUpdateHistory ... update history variables after leaving the routine
    //! @param rConstitutiveInput ... input to the constitutive law (strain, temp gradient etc.)
    //! @param rConstitutiveOutput ... output to the constitutive law (stress, stiffness, heat flux etc.)
    NuTo::Error::eError Evaluate1D(ElementBase* rElement, int rIp,
    		const std::map<NuTo::Constitutive::Input::eInput, const ConstitutiveInputBase*>& rConstitutiveInput,
    		std::map<NuTo::Constitutive::Output::eOutput, ConstitutiveOutputBase*>& rConstitutiveOutput) override;

    //! @brief ... evaluate the constitutive relation in 2D
    //! @param rElement ... element
    //! @param rIp ... integration point
    //! @param rUpdateHistory ... update history variables after leaving the routine
    //! @param rConstitutiveInput ... input to the constitutive law (strain, temp gradient etc.)
    //! @param rConstitutiveOutput ... output to the constitutive law (stress, stiffness, heat flux etc.)
    NuTo::Error::eError Evaluate2D(ElementBase* rElement, int rIp,
    		const std::map<NuTo::Constitutive::Input::eInput, const ConstitutiveInputBase*>& rConstitutiveInput,
    		std::map<NuTo::Constitutive::Output::eOutput, ConstitutiveOutputBase*>& rConstitutiveOutput) override;

    //! @brief ... evaluate the constitutive relation in 3D
    //! @param rElement ... element
    //! @param rIp ... integration point
    //! @param rUpdateHistory ... update history variables after leaving the routine
    //! @param rConstitutiveInput ... input to the constitutive law (strain, temp gradient etc.)
    //! @param rConstitutiveOutput ... output to the constitutive law (stress, stiffness, heat flux etc.)
    NuTo::Error::eError Evaluate3D(ElementBase* rElement, int rIp,
    		const std::map<NuTo::Constitutive::Input::eInput, const ConstitutiveInputBase*>& rConstitutiveInput,
    		std::map<NuTo::Constitutive::Output::eOutput, ConstitutiveOutputBase*>& rConstitutiveOutput) override;

    //! @brief ... create new static data object for an integration point
    //! @return ... pointer to static data object
    ConstitutiveStaticDataBase* AllocateStaticDataEngineeringStress_EngineeringStrain1D(const ElementBase* rElement) const override;

    //! @brief ... create new static data object for an integration point
    //! @return ... pointer to static data object
    ConstitutiveStaticDataBase* AllocateStaticDataEngineeringStress_EngineeringStrain2D(const ElementBase* rElement) const override;

    //! @brief ... create new static data object for an integration point
    //! @return ... pointer to static data object
    ConstitutiveStaticDataBase* AllocateStaticDataEngineeringStress_EngineeringStrain3D(const ElementBase* rElement) const override;

    //! @brief ... performs the return mapping procedure in 3D
    //! @param rElement ... structure
    //! @param rIp ... integration point
    //! @param rEngineeringStrain ... engineering strain
    //! @param rNewStress ... new stress (if a 0-pointer is given, no values are written)
    //! @param rNewTangent ... new tangent matrix (if a 0-pointer is given, no values are written)
    //! @param rNewStaticData ... new static data (if a 0-pointer is given, no values are written)
    NuTo::Error::eError ReturnMapping3D(const ElementBase* rElement,int rIp,
    		const EngineeringStrain3D& rEngineeringStrain,
    		EngineeringStress3D* rNewStress,
    		ConstitutiveTangentLocal<6,6>* rNewTangent,
    		ConstitutiveStaticDataDamageViscoPlasticity3D* rNewStaticData,
    		Logger& rLogger)const;

    ///////////////////////////////////////////////////////////////////////////


    // calculate coefficients of the material matrix
    void CalculateCoefficients2DPlainStress(double& C11, double& C12, double& C33) const;
    void CalculateCoefficients3D(double& C11, double& C12, double& C44) const;

    // parameters /////////////////////////////////////////////////////////////
    //! @brief ... get density
    //! @return ... density
    virtual double GetDensity() const override;

    //! @brief ... set density
    //! @param rRho ... density
    virtual void SetDensity(double rRho) override;

    //! @brief ... get Young's modulus
    //! @return ... Young's modulus
    double GetYoungsModulus() const override;

    //! @brief ... set Young's modulus
    //! @param rE ... Young's modulus
    void SetYoungsModulus(double rE) override;

    //! @brief ... get Poisson's ratio
    //! @return ... Poisson's ratio
    double GetPoissonsRatio() const override;

    //! @brief ... set Poisson's ratio
    //! @param rNu ... Poisson's ratio
    void SetPoissonsRatio(double rNu) override;

    //! @brief ... get thermal expansion coefficient
    //! @return ... thermal expansion coefficient
    double GetThermalExpansionCoefficient() const override;

    //! @brief ... set thermal expansion coefficient
    //! @param rAlpha ... thermal expansion coefficient
    void SetThermalExpansionCoefficient(double rNu) override;

    //! @brief ... get tensile strength
    //! @return ... tensile strength
    double GetTensileStrength() const;

    //! @brief ... set tensile strength
    //! @param rTensileStrength...  tensile strength
    void SetTensileStrength(double rTensileStrength);

    //! @brief ... get compressive strength
    //! @return ... compressive strength
    double GetCompressiveStrength() const;

    //! @brief ... set compressive strength
    //! @param rCompressiveStrength...  compressive strength
    void SetCompressiveStrength(double rCompressiveStrength);

    //! @brief ... get biaxial compressive strength
    //! @return ... biaxial compressive strength
    double GetBiaxialCompressiveStrength() const;

    //! @brief ... set biaxial compressive strength
    //! @param rBiaxialCompressiveStrength...  biaxial compressive strength
    void SetBiaxialCompressiveStrength(double rBiaxialCompressiveStrength);

    ///////////////////////////////////////////////////////////////////////////

    //! @brief ... get type of constitutive relationship
    //! @return ... type of constitutive relationship
    //! @sa eConstitutiveType
    Constitutive::eConstitutiveType GetType() const override;

    //! @brief ... check parameters of the constitutive relationship
    void CheckParameters()const;

    //! @brief ... check compatibility between element type and type of constitutive relationship
    //! @param rElementType ... element type
    //! @return ... <B>true</B> if the element is compatible with the constitutive relationship, <B>false</B> otherwise.
    bool CheckElementCompatibility(Element::eElementType rElementType) const override;

    //! @brief ... print information about the object
    //! @param rVerboseLevel ... verbosity of the information
    //! @param rLogger stream for the output
    void Info(unsigned short rVerboseLevel, Logger& rLogger) const override;

    //! @brief ... returns true, if a material model has tmp static data (which has to be updated before stress or stiffness are calculated)
    //! @return ... see brief explanation
    bool HaveTmpStaticData() const override
    {
    	return false;
    }

#ifdef ENABLE_SERIALIZATION
    //! @brief serializes the class
    //! @param ar         archive
    //! @param version    version
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version);
#endif // ENABLE_SERIALIZATION

protected:
    //! @brief ... Young's modulus \f$ E \f$
    double mE;

    //! @brief ... Poisson's ratio \f$ \nu \f$
    double mNu;

    //! @brief ... density \f$ \rho \f$
    double mRho;

    //! @brief ... thermal expansion coefficient \f$ \alpha \f$
    double mThermalExpansionCoefficient;

    //! @param ... Tensile Strength \f$ \TensileStrength \f$
    double mTensileStrength;

    //! @param ... Compressive Strength \f$ \CompressiveStrength \f$
    double mCompressiveStrength;

    //! @param ... Biaxial Compressive Strength \f$ \BiaxialCompressiveStrength \f$
    double mBiaxialCompressiveStrength;

    //! @brief ... check if density is positive
    //! @param rRho ... density
    void CheckDensity(double rRho) const;

    //! @brief ... check if Young's modulus is positive
    //! @param rE ... Young's modulus
    void CheckYoungsModulus(double rE) const;

    //! @brief ... check if Poisson's ratio is valid
    //! @param rNu ... Poisson's ratio
    void CheckPoissonsRatio(double rNu) const;

    //! @brief ... check thermal expansion coefficient
    //! @param rAlpha ... thermal expansion coefficient
    void CheckThermalExpansionCoefficient(double rAlpha) const;

    //! @brief ... check if tensile strength is positive
    //! @param rTensileStrength ...
    void CheckTensileStrength(double rTensileStrength) const;

    //! @brief ... check if compressive strength is positive
    //! @param rRadius ... compressive strength
    void CheckCompressiveStrength(double rCompressiveStrength) const;

    //! @brief ... check if biaxial compressive strength is positive
    //! @param rBiaxialCompressiveStrength ... biaxial compressive strength
    void CheckBiaxialCompressiveStrength(double rBiaxialCompressiveStrength) const;

private:

    //! @brief ... calculates the residual vector of an incremental formulation
    //! @brief ... rElasticEngineeringStrain ... elastic engineering strain at the end of time increment
    NuTo::FullVector<double,Eigen::Dynamic> Residual(NuTo::FullVector<double,Eigen::Dynamic> rz) const
//    		const EngineeringStrain3D& rElasticEngineeringStrain) const
		{
			NuTo::FullVector<double,Eigen::Dynamic> residual(rz.GetNumRows());

			residual[0] = 2.*rz[1]*rz[2] + 2.*rz[0]*rz[1];
			residual[1] = rz[1]*rz[1] + 2.*rz[0]*rz[2];
			residual[2] = rz[1]*rz[1] + rz[2]*rz[2] -3.;

			return residual;
		}

    //! @brief ... calculates the analytical Jacobi matrix:= derivative of the Residual
    //! @brief ... rElasticEngineeringStrain ... elastic engineering strain at the end of time increment
    NuTo::FullMatrix<double,Eigen::Dynamic,Eigen::Dynamic> DResidualAn(NuTo::FullVector<double,Eigen::Dynamic> rz) const
//    		const EngineeringStrain3D& rElasticEngineeringStrain) const
		{
			NuTo::FullMatrix<double,Eigen::Dynamic,Eigen::Dynamic> deriv(rz.GetNumRows(),rz.GetNumRows());

			deriv(0,0) = 2*rz[1];
			deriv(0,1) = 2*rz[0] + 2*rz[2];
			deriv(0,2) = 2*rz[1];

			deriv(1,0) = 2*rz[2];
			deriv(1,1) = 2*rz[1];
			deriv(1,2) = 2*rz[0];

			deriv(2,0) = 0.;
			deriv(2,1) = 2*rz[1];
			deriv(2,2) = 2*rz[2];

			return deriv;
		}

    //! @brief ... calculates 0.5*Residual^2 and updates fvec = Residual(rz)
    double Fmin(NuTo::FullVector<double,Eigen::Dynamic> rz, NuTo::FullVector<double,Eigen::Dynamic> &fvec) const {
    	fvec = this->Residual(rz);
    	double sum=0;
		sum = fvec.dot(fvec);
		return 0.5*sum;
    }

    //! @brief ... calculates the numerical Jacobi matrix:= derivative of the Residual
    //! @brief ... rElasticEngineeringStrain ... elastic engineering strain at the end of time increment
    NuTo::FullMatrix<double,Eigen::Dynamic,Eigen::Dynamic> DResidualNum(NuTo::FullVector<double,Eigen::Dynamic> rz,
    		NuTo::FullVector<double,Eigen::Dynamic> &fvec) const {
    	const double EPS = 1.0e-8;
    	int n=rz.GetNumRows();
    	NuTo::FullMatrix<double,Eigen::Dynamic,Eigen::Dynamic> deriv(n,n);
    	NuTo::FullVector<double,Eigen::Dynamic> xh=rz;
    	for (int j=0;j<n;j++) {
    		double temp=xh[j];
    		double h=EPS*fabs(temp);
    	   	if (h == 0.0) h=EPS;
    	   		xh[j]=temp+h;
    	   		h=xh[j]-temp;
    	   		NuTo::FullVector<double,Eigen::Dynamic> f=this->Residual(xh);
    	   		xh[j]=temp;
    	   		for (int i=0;i<n;i++)
    	   			deriv(i,j)=(f[i]-fvec[i])/h;
    	   		}
    	return deriv;
    }

    //! @brief ... the routine performs line search correction of the Newton step
// NR    template <class T>
    void LineSearch(NuTo::FullVector<double,Eigen::Dynamic> &xold,
    		const double fold,
    		NuTo::FullVector<double,Eigen::Dynamic> &g,
    		NuTo::FullVector<double,Eigen::Dynamic> &p,
    		NuTo::FullVector<double,Eigen::Dynamic> &x,
// NR    		double &f, const double stpmax, bool &check, T &func) {
		double &f, const double stpmax, bool &check, NuTo::FullVector<double,Eigen::Dynamic> &fvec) {  // AnstattNR
    	const double ALF=1.0e-4, TOLX=numeric_limits<double>::epsilon();
    	double a,alam,alam2=0.0,alamin,b,disc,f2=0.0;
    	double rhs1,rhs2,slope=0.0,sum=0.0,test,tmplam;

    	check = false;
    //                                              // OPTIMIZED
    //	for (i=0;i<n;i++) sum += p[i]*p[i];   		// OPTIMIZED
    //	sum=sqrt(sum);                        		// OPTIMIZED
    	sum = p.lpNorm<2>();                  		// OPTIMIZED

    	if (sum > stpmax)
    //  											// OPTIMIZED
    //		for (i=0;i<n;i++)                       // OPTIMIZED
    //			p[i] *= stpmax/sum;                 // OPTIMIZED
    		p *= stpmax/sum;                        // OPTIMIZED

    //                                              // OPTIMIZED
    //	for (i=0;i<n;i++)                           // OPTIMIZED
    //		slope += g[i]*p[i];                     // OPTIMIZED
    	slope = g.dot(p);	                        // OPTIMIZED

    	if (slope >= 0.0){
    		cout << "Roundoff problem in lnsrch." << endl;
    		throw("Roundoff problem in lnsrch.");
    	}
    	test=0.0;
    //                 								// OPTIMIZED
    //	for (i=0;i<n;i++) {
    //		temp=fabs(p[i])/max(fabs(xold[i]),1.0); // OPTIMIZED
    //		if (temp > test) test=temp;             // OPTIMIZED
    //	}    										// OPTIMIZED
    	test = ( p.array().abs() / xold.array().abs().max(1.0) ).maxCoeff();  // OPTIMIZED

    	alamin=TOLX/test;
    	alam=1.0;
    	for (;;) {
    // 												// OPTIMIZED
    //		for (i=0;i<n;i++) x[i]=xold[i]+alam*p[i];  // OPTIMIZED
    		x = xold + alam*p; 						// OPTIMIZED
    // NR		f=func(x);
    		f = this->Fmin(x,fvec);
    		if (alam < alamin) {
    //												// OPTIMIZED
    //			for (i=0;i<n;i++) x[i]=xold[i]; 	// OPTIMIZED
    			x = xold;							// OPTIMIZED
    			check=true;
    			return;
    		} else if (f <= fold+ALF*alam*slope) return;
    		else {
    			if (alam == 1.0)
    				tmplam = -slope/(2.0*(f-fold-slope));
    			else {
    				rhs1=f-fold-alam*slope;
    				rhs2=f2-fold-alam2*slope;
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
    		f2 = f;
    		alam=max(tmplam,0.1*alam);
    	}
    }

    //! @brief ... the routine performs Newton-Raphson integration
//NR    template <class T>
    void Newton(NuTo::FullVector<double,Eigen::Dynamic> &x, bool &check,
//NR    		T &vecfunc,
   		NuTo::FullMatrix<double,Eigen::Dynamic,Eigen::Dynamic>
    	(DamageViscoPlasticityEngineeringStress::*fdjacAn)(NuTo::FullVector<double,Eigen::Dynamic>) const = 0) {
    	const int MAXITS=200;
    	const double TOLF=1.0e-8,TOLMIN=1.0e-12,STPMX=100.0;
    	const double TOLX=numeric_limits<double>::epsilon();
    	int its,n=x.GetNumRows();
    	double den,f,fold,stpmax,test;
    	NuTo::FullVector<double,Eigen::Dynamic> g(n),p(n),xold(n);
    	NuTo::FullMatrix<double,Eigen::Dynamic,Eigen::Dynamic> fjac(n,n);
//NR    	NRfmin<T> fmin(vecfunc);
//NR    	NRfdjac<T> fdjac(vecfunc);
//NR    	NuTo::FullVector<double,Eigen::Dynamic> &fvec=fmin.fvec;
    	NuTo::FullVector<double,Eigen::Dynamic> fvec;   // AnstattNR

//NR    	f=fmin(x);
    	f = this->Fmin(x, fvec); // AnstattNR
    //												// OPTIMIZED
    //	test=0.0;
    //	for (i=0;i<n;i++)							// OPTIMIZED
    //		if (fabs(fvec[i]) > test) test=fabs(fvec[i]);   // OPTIMIZED
    	test = fvec.array().abs().maxCoeff();		// OPTIMIZED
    	if (test < 0.01*TOLF) {
    		check=false;
    		return;
    	}
    //												// OPTIMIZED
    //	sum=0.0;									// OPTIMIZED
    //	for (i=0;i<n;i++) sum += x[i]*x[i];			// OPTIMIZED
    //	stpmax=STPMX*max(sqrt(sum),double(n));		// OPTIMIZED
    	stpmax=STPMX*max(x.lpNorm<2>(),double(n));	// OPTIMIZED

    	for (its=0;its<MAXITS;its++) {
    		cout<< "===== Iteration =====" << its <<endl;   // Test
    		if (fdjacAn != 0) {         			// if analytical Jacobi is there
    			fjac = (this->*fdjacAn)(x);   			// analytical Jacobi matrix
    			cout<<"*** Analytical ***"<<endl;  	// Test
    			cout << fjac << endl;             	// Test
    		} else {                   				// if not analytic -> take numeric
    			fjac=this->DResidualNum(x,fvec);    // numerical Jacobi matrix
    			cout<<"*** Numerical ***"<<endl;   	// Test
    			cout << fjac << endl;	           	// Test
    		}
    //												// OPTIMIZED
    //		for (i=0;i<n;i++) {						// OPTIMIZED
    //			sum=0.0;							// OPTIMIZED
    //			for (j=0;j<n;j++) sum += fjac(j,i)*fvec[j];   // OPTIMIZED
    //			g[i]=sum;							// OPTIMIZED
    //		}										// OPTIMIZED
    		g = fjac.transpose() * fvec;			// OPTIMIZED
    //												// OPTIMIZED
    //		for (i=0;i<n;i++) xold[i]=x[i];         // OPTIMIZED
    		xold = x;
    		fold=f;
    // 												// OPTIMIZED
    //		for (i=0;i<n;i++) p[i] = -fvec[i];      // OPTIMIZED
    		p = -fvec;								// OPTIMIZED
    //												// Test
    		cout << "fvec = " << fvec.transpose()<<endl;   // Test
    		cout << "p = " << fjac.fullPivLu().solve(-fvec).transpose()<<endl;   // Test

    //												// OPTIMIZED
    //		LUdcmp alu(fjac);						// OPTIMIZED
    //		alu.solve(p,p)							// OPTIMIZED
    		p =	fjac.fullPivLu().solve(-fvec).transpose();		// LU SOLVER of fjac * p = -fvec
    //															// SVD SOLVER
    //		p = fjac.jacobiSvd().solve(-fvec).transpose();      // SVD SOLVER
//NR    		LineSearch(xold,fold,g,p,x,f,stpmax,check,fmin);
    		LineSearch(xold,fold,g,p,x,f,stpmax,check,fvec);    // AnstattNR

    //												// OPTIMIZED
    //		test=0.0;								// OPTIMIZED
    //		for (i=0;i<n;i++)						// OPTIMIZED
    //			if (fabs(fvec[i]) > test) test=fabs(fvec[i]);   // OPTIMIZED
    		test = fvec.array().abs().maxCoeff();	// OPTIMIZED

    		if (test < TOLF) {
    			check=false;
    			return;
    		}

    		if (check) {
    			den=max(f,0.5*n);
    //												// OPTIMIZED
    //			test=0.0;							// OPTIMIZED
    //			for (i=0;i<n;i++) {					// OPTIMIZED
    //				temp=fabs(g[i])*max(fabs(x[i]),1.0)/den;   // OPTIMIZED
    //				if (temp > test) test=temp;		// OPTIMIZED
    //			}									// OPTIMIZED
    			test = ( g.array().abs() * x.array().abs().max(1.0) ).maxCoeff() / den;   // OPTIMIZED
    			check=(test < TOLMIN);
    			return;
    		}
    //												// OPTIMIZED
    //		test=0.0;								// OPTIMIZED
    //		for (i=0;i<n;i++) {						// OPTIMIZED
    //			temp=(fabs(x[i]-xold[i]))/max(fabs(x[i]),1.0);   // OPTIMIZED
    //			if (temp > test) test=temp;			// OPTIMIZED
    //		}										// OPTIMIZED
    		test = ( (x.array()-xold.array()).abs() / x.array().abs().max(1.0) ).maxCoeff();   // OPTIMIZED

    		if (test < TOLX)
    			return;
    	}
    	throw("MAXITS exceeded in newt");
    }

};


}


#ifdef ENABLE_SERIALIZATION
BOOST_CLASS_EXPORT_KEY(NuTo::DamageViscoPlasticityEngineeringStress)
#endif //ENABLE_SERIALIZATION

#endif // DamageViscoPlasticityEngineeringStress_H_
