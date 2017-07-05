
#pragma once

#include "mechanics/timeIntegration/NewmarkBase.h"
#include "mechanics/timeIntegration/NewmarkDirect.h"
#include <boost/math/constants/constants.hpp>
#include <math.h>

namespace NuTo
{
//! @author Vitaliy M. Kindrachuk
//! @date December 2014
//! @brief ... standard class for implicit time integration including cycle jump based on a Fourier formulation
class JumpDirect : public NewmarkDirect
{

public:
    //! @brief constructor
    JumpDirect(StructureBase* rStructure);

    void SetMinLineSearchStep(double rMinLineSearchStep)
    {
        mMinLineSearchStep = rMinLineSearchStep;
    }

    double GetMinLineSearchStep() const
    {
        return mMinLineSearchStep;
    }

    //!@brief sets the number of increments of a single cycle, this should be divisible by 4, because the cycle consist
    //!of 4 symmetric parts
    //!@param sets the number of increments of a single cycle, this should be divisible by 4, because the cycle consist
    //!of 4 symmetric parts
    void SetHarmonicIncrementation(int rHarmonicIncrementation)
    {
        if (fmod(rHarmonicIncrementation, 4) == 0)
        {
            mHarmonicIncrementation = rHarmonicIncrementation;
        }
        else
        {
            mHarmonicIncrementation = (rHarmonicIncrementation - fmod(rHarmonicIncrementation, 4) + 4);
        }
    }

    //!@brief gets the number of increments of a single cycle, this should be divisible by 4, because the cycle consist
    //!of 4 parts
    //!@param gets the number of increments of a single cycle, this should be divisible by 4, because the cycle consist
    //!of 4 parts
    int GetHarmonicIncrementation() const
    {
        return mHarmonicIncrementation;
    }

    //! @brief if mHarmonicExtrapolation == false then performs a Newmark, otherwise performs Newmark till 3 cycles and
    //! then cycle jump
    //! @param if mHarmonicExtrapolation == false then performs a Newmark, otherwise performs Newmark till 3 cycles and
    //! then cycle jump
    void SetHarmonicExtrapolation(bool rHarmonicExtrapolation, double* rHarmonicExtrapolationTolerance = 0)
    {
    	if (mHarmonicExcitation == true) {
    		mHarmonicExtrapolation = rHarmonicExtrapolation;
    		if (mHarmonicExtrapolation && rHarmonicExtrapolationTolerance != 0) {
    			mHarmonicExtrapolationTolerance = *rHarmonicExtrapolationTolerance;
			}
		} else {
			throw Exception("[NuTo::JumpDirect::SetHarmonicConstraint] Set HarmonicConstraint at first!");
		}

    }

    //! @brief if mHarmonicExtrapolation == false then performs a Newmark, otherwise performs Newmark till 3 cycles and
    //! then cycle jump
    //! @param if mHarmonicExtrapolation == false then performs a Newmark, otherwise performs Newmark till 3 cycles and
    //! then cycle jump
    bool GetHarmonicExtrapolation() const
    {
        return mHarmonicExtrapolation;
    }

    //!@brief sets the harmonic excitations in the format {DisplacementAmplitude,Frequency,NumberOfCycles} if
    //!displacement controlled or
    //!@brief in the format {LoadingAmplitude,Frequency,Number OfCycles} if load controlled
    //!@param sets the harmonic excitations in the format {DisplacementAmplitude,Frequency,NumberOfCycles} or
    //!{LoadingAmplitude,Frequency,Number OfCycles}
    // rHarmonicConstraintFactor[*,0] = displacement or load amplitude, arbitrary sign
    // rHarmonicConstraintFactor[*,1] = frequency, positive
    // rHarmonicConstraintFactor[*,2] = number of cycles, at least 1
    // The harmonic excitations are applied to the mTimeDependentConstraint !
    // Currently is only a monoharmonic constraint applied, that is rHarmonicFactor.rows()=1
    void SetHarmonicExcitation(const NuTo::FullMatrix<double, Eigen::Dynamic, 3>& rHarmonicFactor)
    {
    	if (mTimeDependentConstraint == -1 && mTimeDependentLoadCase == -1) {
    		throw Exception("[NuTo::JumpDirect::SetHarmonicExcitation] the harmonic excitation can be only applied to the time dependent constraint or load case.");
		}

    	if (mTimeDependentConstraint != -1 && mTimeDependentLoadCase != -1) {
    		throw Exception("[NuTo::JumpDirect::SetHarmonicExcitation] the harmonic excitation is currently implemented for either time dependent constraint or load case.");
		}

    	if (rHarmonicFactor.rows()>1)
    		throw Exception("[NuTo::JumpDirect::SetHarmonicExcitation] currently a monoharmonic constraint is implemented, number of rows must be 1.");

    	//check, whether frequencies and number of cycles are positive
    	for (int count=0; count<rHarmonicFactor.rows(); count++)
    	{
    		if (rHarmonicFactor(count,1)<=0.)
    			throw Exception("[NuTo::JumpDirect::SetHarmonicExcitation] the frequency should always be positive.");
    		if (rHarmonicFactor(count,2)<=1.)
    			throw Exception("[NuTo::JumpDirect::SetHarmonicExcitation] number of cycles should be at least 1.");
    	}

    	mHarmonicExcitation = true;
    	mHarmonicFactor = rHarmonicFactor;
    }

    //! @brief returns true, if there is a harmonic excitation applied
    bool GetHarmonicExcitation() const
    {
        return mHarmonicExcitation;
    }

    //! @brief apply the new rhs of the constraints as a function of the current time delta
    //! @brief the time dependent constraint is followed by a sine excitation
    double CalculateTimeDependentConstraintFactor(double curTime)
    {
    	if (mTimeDependentConstraintFactor.rows()==0)
    		throw Exception("[NuTo::JumpDirect::CalculateTimeDependentConstraintFactor] the harmonic excitation can be only applied to the time dependent constraint.");

    	// calculate the end time of the time dependent constraint. After this time the constraint is harmonic
    	double timeDependentConstraintTime(mTimeDependentConstraintFactor(mTimeDependentConstraintFactor.rows()-1,0));
		double timeDependentConstraintFactor;

    	if (curTime <= timeDependentConstraintTime) {
        	// calculate timeDependentConstraintFactor during the linear time dependent constraint
    		return NewmarkDirect::CalculateTimeDependentConstraintFactor(curTime);
		} else {
			// time dependent constraint is continued by the harmonic excitation

			// calculate constraint at the end time of the time dependent constraint
			timeDependentConstraintFactor = NewmarkDirect::CalculateTimeDependentConstraintFactor(timeDependentConstraintTime);

			if (mHarmonicFactor.rows()==0) {
				return timeDependentConstraintFactor;
			}

			// add harmonic excitation
			double frequency(mHarmonicFactor(0,1));
			double amplitude(mHarmonicFactor(0,0));
			const double pi = boost::math::constants::pi<double>();
			timeDependentConstraintFactor += amplitude*sin(2*pi*frequency*(curTime - timeDependentConstraintTime));
			return timeDependentConstraintFactor;
		}
    }

    //! @brief calculates the applied load cases in a single load vector as a function of the current time delta
    //! @brief the time dependent linear load is followed by a sine excitation
    void CalculateExternalLoad(StructureBase& rStructure, double curTime, Eigen::VectorXd& rLoad_j,
                               Eigen::VectorXd& rLoad_k)
    {
    	rLoad_j.Resize(rStructure.GetNumActiveDofs());
    	rLoad_k.Resize(rStructure.GetNumDofs()-rStructure.GetNumActiveDofs());

    	if (mTimeDependentLoadCase == -1) {
			// there is any time dependent load case
    		NewmarkDirect::CalculateExternalLoad(rStructure,curTime,rLoad_j,rLoad_k);
    	} else {
			// there is a time dependent load case defined

			// check that the TimeDependentLoadFactor is set
			if (mTimeDependentLoadFactor.rows()==0)
			{
				throw Exception("[NuTo::JumpDirect::CalculateExternalLoad] TimeDependentLoadFactor not set.");
			}
			// extract excitation parameters
			double frequency(mHarmonicFactor(0,1));
			double amplitude(mHarmonicFactor(0,0));
			const double pi = boost::math::constants::pi<double>();

			// calculate load amplitude
			Eigen::VectorXd LoadAmplitude_j,LoadAmplitude_k;
			LoadAmplitude_j.Resize(rStructure.GetNumActiveDofs());
			LoadAmplitude_k.Resize(rStructure.GetNumDofs()-rStructure.GetNumActiveDofs());

			LoadAmplitude_j.setOnes(); LoadAmplitude_k.setOnes();   // check this function in inverseMatrixExample
			LoadAmplitude_j *= amplitude;
			LoadAmplitude_k *= amplitude;

	    	// calculate the end time of the time dependent load case. After this time the load is harmonic
	    	double timeDependentLoadCaseTime(mTimeDependentLoadFactor(mTimeDependentLoadFactor.rows()-1,0));

	    	if (curTime <= timeDependentLoadCaseTime) {
	        	// calculate timeDependentLoadFactor during the linear time dependent load
	    		NewmarkDirect::CalculateExternalLoad(rStructure,curTime,rLoad_j,rLoad_k);
	    	} else {
				// time dependent load case is continued by the harmonic excitation
				double s(mTimeDependentLoadFactor(mTimeDependentLoadFactor.rows()-1,1));	// maximal linear load factor

//				rLoad_j=mLoadVectorStatic_j+mLoadVectorTimeDependent_j*s + LoadAmplitude_j*sin(2*pi*frequency*(curTime - timeDependentLoadCaseTime));
//				rLoad_k=mLoadVectorStatic_k+mLoadVectorTimeDependent_k*s + LoadAmplitude_k*sin(2*pi*frequency*(curTime - timeDependentLoadCaseTime));
				rLoad_j=mLoadVectorStatic_j+mLoadVectorTimeDependent_j*s + mLoadVectorTimeDependent_j*amplitude*sin(2*pi*frequency*(curTime - timeDependentLoadCaseTime));
				rLoad_k=mLoadVectorStatic_k+mLoadVectorTimeDependent_k*s + mLoadVectorTimeDependent_k*amplitude*sin(2*pi*frequency*(curTime - timeDependentLoadCaseTime));
	    	}
		}
    }

    void SetTimeAndTimeStep(double& curTime, double& timeStep, double rTimeDelta)
    {
        // calculate time and time step close to the end time of the integration interval
        if (rTimeDelta - curTime < 0.2 * timeStep)
        {
            timeStep += rTimeDelta - curTime;
            curTime = rTimeDelta;
            return;
        }

        // calculate time and time step during the time dependent constraint (where the linear interpolation is
        // performed)
        if (mTimeDependentConstraintFactor.rows() != 0 && mAutomaticTimeStepping)
        {
            int curStep(0);
            while (mTimeDependentConstraintFactor(curStep, 0) < curTime &&
                   curStep < mTimeDependentConstraintFactor.rows() - 1)
                curStep++;
            if (curStep == 0)
                curStep++;

            // extract the two data points
            double t1 = mTimeDependentConstraintFactor(curStep - 1, 0);
            double t2 = mTimeDependentConstraintFactor(curStep, 0);

            // if curTime is close to the end of the time increment then set it to the end
            if (t2 - curTime < 0.5 * timeStep && curTime < t2)
            {
                timeStep += t2 - curTime;
                curTime = t2;
                return;
            }
            // if curTime jumps over a data point, then set it to this data point
            if (curTime - timeStep < t1 - 0.2 * timeStep)
            {
                timeStep -= curTime - t1;
                curTime = t1;
                return;
            }
            // if curTime jumps over the end of the time dependent constraint, then set it
            // to this end (this point is at the same time the beginning of the harmonic excitation)
            if (curTime - timeStep < t2 - 0.2 * timeStep && curTime > t2)
            {
                timeStep -= curTime - t2;
                curTime = t2;
                return;
            }
            if (curTime < t2)
            {
                return;
            }
        }

        // calculate time and time step during the time dependent constraint (where the linear interpolation is
        // performed)
        if (mTimeDependentLoadFactor.rows() != 0 && mAutomaticTimeStepping)
        {
            int curStep(0);
            while (mTimeDependentLoadFactor(curStep, 0) < curTime && curStep < mTimeDependentLoadFactor.rows() - 1)
                curStep++;
            if (curStep == 0)
                curStep++;

            // extract the two data points
            double t1 = mTimeDependentLoadFactor(curStep - 1, 0);
            double t2 = mTimeDependentLoadFactor(curStep, 0);

            // if curTime is close to the end of the time increment then set it to the end
            if (t2 - curTime < 0.5 * timeStep && curTime < t2)
            {
                timeStep += t2 - curTime;
                curTime = t2;
                return;
            }
            // if curTime jumps over a data point, then set it to this data point
            if (curTime - timeStep < t1 - 0.2 * timeStep)
            {
                timeStep -= curTime - t1;
                curTime = t1;
                return;
            }
            // if curTime jumps over the end of the time dependent constraint, then set it
            // to this end (this point is at the same time the beginning of the harmonic excitation)
            if (curTime - timeStep < t2 - 0.2 * timeStep && curTime > t2)
            {
                timeStep -= curTime - t2;
                curTime = t2;
                return;
            }
            if (curTime < t2)
            {
                return;
            }
        }

        // calculate time step during harmonic excitation, the number of increments should be divisible by 4
        if (mHarmonicExcitation == true &&
            std::abs(timeStep - 1. / (mHarmonicIncrementation * mHarmonicFactor(0, 1))) > 0.001 * timeStep)
        {
            curTime += 1. / (mHarmonicIncrementation * mHarmonicFactor(0, 1)) - timeStep;
            timeStep = 1. / (mHarmonicIncrementation * mHarmonicFactor(0, 1));
            return;
        }
    }

    //! @brief returns true, if the method is only conditionally stable (for unconditional stable, this is false)
    bool HasCriticalTimeStep() const
    {
        return false;
    }

    //! @brief calculate the critical time step for explicit routines
    //! for implicit routines, this will simply return zero (cmp HasCriticalTimeStep())
    double CalculateCriticalTimeStep() const
    {
        return 0;
    }
    //!@brief calculate the global stiffness matrix for the BVP with mean displacement (rFourier = 0) or with
    //!displacement amplitude (rFourier = 1)
    void CalculateGlobalModifiedStiffness(NuTo::SparseMatrixCSRVector2General<double>* rStiffnessMod, int rFourierMode);

    //!@brief find equilibrium by calculating the Fourier coefficients (which are the displacement fields)
    // rDisp_Mean_j ... mean displacement (active DOF, rFourier = 0)
    // rDisp_Mean_k ... mean displacement (dependent DOF, rFourier = 0)
    // rDisp_Ampl_j ... displacement amplitude (active DOF, rFourier = 1)
    // rDisp_Mean_k ... displacement amplitude (dependent DOF, rFourier = 1)
    void CalculateFourierCoefficients(Eigen::VectorXd* rDisp_Mean_j, Eigen::VectorXd* rDisp_Mean_k,
                                      Eigen::VectorXd* rDisp_Ampl_j, Eigen::VectorXd* rDisp_Ampl_k,
                                      const Eigen::VectorXd* rIntForce_Mean_j, const Eigen::VectorXd* rIntForce_Mean_k,
                                      const Eigen::VectorXd* rIntForce_Max_j, const Eigen::VectorXd* rIntForce_Max_k);

    //!@brief straight-forward integration of a single cycle with a prescribed Fourier coefficients
    // the displacement fields have the same meaning as above, see CalculateFourierCoefficients
    // rIncludePostProcess ...false, if no postprocessing should be done during integration
    // rIncludePostProcess ...true, postprocessing will be done
    // Postprocessing should be done if the jump is acceptable; in this case the IntegrateSinleCycle should be repeated
    // with a true option
    void IntegrateSingleCycle(Eigen::VectorXd* rDisp_Mean_j, Eigen::VectorXd* rDisp_Mean_k,
                              Eigen::VectorXd* rDisp_Ampl_j, Eigen::VectorXd* rDisp_Ampl_k, bool rIncludePostProcess);

    //!@brief updates DofTypes after the CalculateFourierCoefficients routine
    // The CalculateFourierCoefficients routine determines the DISPLACEMENTS Dof only. For a coupled problem
    // DISPLACEMENTS/DofType, the another DofType
    // has to be updated too. This routine updates the DofTypesfor the following list of Dofs:
    // NONLOCALEQSTRAIN
    // rDisp_Mean_j ... mean displacement (active DOF, rFourier = 0)
    // rDisp_Mean_k ... mean displacement (dependent DOF, rFourier = 0)
    // rDisp_Ampl_j ... displacement amplitude (active DOF, rFourier = 1)
    // rDisp_Mean_k ... displacement amplitude (dependent DOF, rFourier = 1)
    void CalculateFourierCoefficientsCoupledDofs(Eigen::VectorXd* rDisp_Mean_j, Eigen::VectorXd* rDisp_Mean_k,
                                                 Eigen::VectorXd* rDisp_Ampl_j, Eigen::VectorXd* rDisp_Ampl_k);

    //! @brief performs the time integration
    //! @param rTimeDelta ... length of the simulation
    NuTo::Error::eError Solve(double rTimeDelta);

    //! @brief ... Info routine that prints general information about the object (detail according to verbose level)
    void Info() const;

    //! @brief ... Return the name of the class, this is important for the serialize routines, since this is stored in
    //! the file
    //!            in case of restoring from a file with the wrong object type, the file id is printed
    //! @return    class name
    virtual std::string GetTypeId() const;

protected:
    double mMinLineSearchStep;
    bool mHarmonicExcitation;
    bool mHarmonicExtrapolation;
    double mHarmonicExtrapolationTolerance;
    int mHarmonicIncrementation;
    NuTo::FullMatrix<double, Eigen::Dynamic, 3> mHarmonicFactor;
};
} // namespace NuTo
