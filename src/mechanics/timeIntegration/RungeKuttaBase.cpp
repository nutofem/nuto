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

#include "math/SparseMatrixCSRVector2.h"
#include "mechanics/nodes/NodeBase.h"
#include "mechanics/nodes/NodeEnum.h"
#include "mechanics/groups/Group.h"
#include "mechanics/structures/StructureBase.h"
#include "mechanics/structures/StructureOutputBlockMatrix.h"
#include "mechanics/timeIntegration/RungeKuttaBase.h"
#include "mechanics/timeIntegration/TimeIntegrationEnum.h"

#include "base/Timer.h"

//! @brief constructor
//! @param mDimension number of nodes
NuTo::RungeKuttaBase::RungeKuttaBase (StructureBase* rStructure)  : TimeIntegrationBase (rStructure)
{
    mTimeStep = 0.;
}


//! @brief ... Info routine that prints general information about the object (detail according to verbose level)
void NuTo::RungeKuttaBase::Info()const
{
	TimeIntegrationBase::Info();
}

#ifdef ENABLE_SERIALIZATION
// serializes the class
template void NuTo::RungeKuttaBase::serialize(boost::archive::binary_oarchive & ar, const unsigned int version);
template void NuTo::RungeKuttaBase::serialize(boost::archive::xml_oarchive & ar, const unsigned int version);
template void NuTo::RungeKuttaBase::serialize(boost::archive::text_oarchive & ar, const unsigned int version);
template void NuTo::RungeKuttaBase::serialize(boost::archive::binary_iarchive & ar, const unsigned int version);
template void NuTo::RungeKuttaBase::serialize(boost::archive::xml_iarchive & ar, const unsigned int version);
template void NuTo::RungeKuttaBase::serialize(boost::archive::text_iarchive & ar, const unsigned int version);
template<class Archive>
void NuTo::RungeKuttaBase::serialize(Archive & ar, const unsigned int version)
{
    #ifdef DEBUG_SERIALIZATION
        std::cout << "start serialization of RungeKuttaBase" << "\n";
    #endif
        ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(TimeIntegrationBase)
           & BOOST_SERIALIZATION_NVP(mTimeStep);
    #ifdef DEBUG_SERIALIZATION
        std::cout << "finish serialization of RungeKuttaBase" << "\n";
    #endif
}

#endif // ENABLE_SERIALIZATION


//! @brief perform the time integration
//! @param rStructure ... structure
//! @param rTimeDelta ... length of the simulation
void NuTo::RungeKuttaBase::Solve(double rTimeDelta)
{
    NuTo::Timer timer(__FUNCTION__, mStructure->GetShowTime(), mStructure->GetLogger());

    mStructure->NodeBuildGlobalDofs(__PRETTY_FUNCTION__);

    if (mStructure->GetDofStatus().HasInteractingConstraints())
        throw MechanicsException(__PRETTY_FUNCTION__, "not implemented for constrained systems including multiple dofs.");

    if (mTimeStep==0.)
    {
        if (this->HasCriticalTimeStep())
        {
            mTimeStep = this->CalculateCriticalTimeStep();
        }
        else
        {
            throw MechanicsException("[NuTo::RungeKuttaBase::Solve] time step not set for unconditional stable algorithm.");
        }
    }

    std::cout << "modify computation of critical time step to include the dependence on the time integration scheme." <<std::endl;
    //calculate instead the smallest eigenfrequency, depending on the time integration this gives the critical time step

    std::cout << "time step " << mTimeStep << std::endl;
    std::cout << "number of time steps " << rTimeDelta/mTimeStep << std::endl;



    StructureOutputBlockVector outOfBalance(mStructure->GetDofStatus(), true);

    CalculateStaticAndTimeDependentExternalLoad();

    //store last converged displacements, velocities and accelerations
    auto dof_dt0 = mStructure->NodeExtractDofValues(0);
    auto dof_dt1 = mStructure->NodeExtractDofValues(1);

    StructureOutputBlockMatrix hessian2 = mStructure->BuildGlobalHessian2Lumped();

    //invert the mass matrix
    hessian2.CwiseInvert();

    double curTime  = 0;
    auto extLoad = CalculateCurrentExternalLoad(curTime);
    auto intForce = mStructure->BuildGlobalInternalGradient();

    std::vector<StructureOutputBlockVector> d_dof_dt0_tmp(this->GetNumStages(), mStructure->GetDofStatus());
    std::vector<StructureOutputBlockVector> d_dof_dt1_tmp(this->GetNumStages(), mStructure->GetDofStatus());

    std::vector<double> stageDerivativeFactor(this->GetNumStages()-1);
    while (curTime < rTimeDelta)
    {
        //calculate for delta_t = 0
        auto dof_dt0_new = dof_dt0;
        auto dof_dt1_new = dof_dt1;

        double prevTime(mTime);
        double prevCurTime(curTime);
        for (int countStage=0; countStage<this->GetNumStages(); countStage++)
        {
            //std::cout << "\n stage weight " << GetStageWeights(countStage) << std::endl;
            double deltaTimeStage = this->GetStageTimeFactor(countStage)*mTimeStep;
            this->GetStageDerivativeFactor(stageDerivativeFactor, countStage);
            auto dof_dt0_tmp = dof_dt0;
            auto dof_dt1_tmp = dof_dt1;
            for (int countStage2=0; countStage2<countStage; countStage2++)
            {
                if (stageDerivativeFactor[countStage2]!=0.)
                {
                    dof_dt0_tmp += d_dof_dt0_tmp[countStage2]*(stageDerivativeFactor[countStage2]);
                    dof_dt1_tmp += d_dof_dt1_tmp[countStage2]*(stageDerivativeFactor[countStage2]);
                }
            }

            if (this->HasTimeChanged(countStage)==true)
            {
                curTime=prevCurTime+deltaTimeStage;
                mTime=prevTime+deltaTimeStage;

                UpdateConstraints(mTime);
               
                extLoad = CalculateCurrentExternalLoad(curTime);
            }
            dof_dt0_tmp.K = mStructure->NodeCalculateDependentDofValues(dof_dt0_tmp.J);
            mStructure->NodeMergeDofValues(0,dof_dt0_tmp);
            mStructure->ElementTotalUpdateTmpStaticData();

            //calculate internal force (with update of history variables = true)
            intForce = mStructure->BuildGlobalInternalGradient();

            //update derivatives (ydot or k1,k2,k3,k4) for Runge Kutta
            d_dof_dt0_tmp[countStage] = dof_dt1_tmp*mTimeStep;
            //std::cout << "d_disp_j_tmp " << d_disp_j_tmp[countStage](0) << std::endl;
            d_dof_dt1_tmp[countStage]  = hessian2*(extLoad-intForce)*mTimeStep;
            //std::cout << "d_vel_j_tmp " << d_vel_j_tmp[countStage](0) << std::endl;
            //std::cout << "norm of acc " << (d_vel_j_tmp).norm() << std::endl;

            dof_dt0_new += d_dof_dt0_tmp[countStage]*(GetStageWeights(countStage));
            //std::cout << "disp_j_new " << disp_j_new(0) << std::endl;
            dof_dt1_new  += d_dof_dt1_tmp[countStage]*(GetStageWeights(countStage));
            //std::cout << "vel_j_new " << vel_j_new(0) << std::endl;
        }

        mTime = prevTime + mTimeStep;
        curTime = prevCurTime + mTimeStep;

        //std::cout << "final disp_j_new " << disp_j_new(0) << std::endl;
        dof_dt0_new.K = mStructure->NodeCalculateDependentDofValues(dof_dt0_new.J);
        mStructure->NodeMergeDofValues(0,dof_dt0_new);
        mStructure->ElementTotalUpdateTmpStaticData();
        mStructure->ElementTotalUpdateStaticData();
        //std::cout << "delta disp between time steps" <<  (disp_j-disp_j_new).norm() << std::endl;
        dof_dt0 = dof_dt0_new;
        dof_dt1 = dof_dt1_new;

        //**********************************************
        //PostProcessing
        //**********************************************
        // outOfBalance_j is automatically zero
        //outOfBalance_j.Resize(intForce_j.rows());
        //the acceleration of the dofs k is given by the acceleration of the rhs of the constraint equation
        //this is calculated using finite differencs
        //make sure to recalculate the internal force and external force (if time factor is not 1)

        //acc_k = (bRHSprev-bRHShalf*2+bRHSend)*(4./(timeStep*timeStep))
        //outOfBalance_k = intForce_k - extForce_k + massMatrix_k.asDiagonal()*acc_k;

        // postprocess data for plotting
        this->PostProcess(extLoad-intForce);
    }
}


#ifdef ENABLE_SERIALIZATION
#ifndef SWIG
BOOST_CLASS_EXPORT_IMPLEMENT(NuTo::RungeKuttaBase)
#endif // SWIG
#endif // ENABLE_SERIALIZATION
