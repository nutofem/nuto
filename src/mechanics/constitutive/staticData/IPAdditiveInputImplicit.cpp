#include "math/SparseDirectSolverMUMPS.h"
#include "math/SparseMatrixCSRGeneral.h"
#include "mechanics/constitutive/staticData/IPAdditiveInputImplicit.h"
#include "mechanics/constitutive/laws/AdditiveInputImplicit.h"
#include "mechanics/constitutive/inputoutput/EngineeringStrain.h"
#include "mechanics/nodes/NodeEnum.h"
#include "mechanics/constitutive/ConstitutiveEnum.h"
using namespace NuTo::Constitutive;


IPAdditiveInputImplicit::IPAdditiveInputImplicit(AdditiveInputImplicit& rLaw, const Data& rData)
    : mLaw(rLaw),
      mData(rData)
{
    mData.AllocateAdditionalData(1);
    for(auto& sublaw : mLaw.mSublaws)
    {
        mSublawIPs.push_back(sublaw->CreateIPLaw().release());
    }
}


IPAdditiveInputImplicit::IPAdditiveInputImplicit(const IPAdditiveInputImplicit& rOther)
    : mLaw(rOther.mLaw),
      mData(rOther.mData),
      mSublawIPs(rOther.mSublawIPs)
{}


std::unique_ptr<IPConstitutiveLawBase> IPAdditiveInputImplicit::Clone() const
{
    return std::make_unique<IPAdditiveInputImplicit>(*this);
}


NuTo::ConstitutiveBase& IPAdditiveInputImplicit::GetConstitutiveLaw() const
{
    return mLaw;
}


void IPAdditiveInputImplicit::AllocateAdditional(int rNum)
{
    for(auto& sublaw : mSublawIPs)
    {
        sublaw.AllocateAdditional(rNum);
    }
}


void IPAdditiveInputImplicit::ShiftToPast()
{
    for(auto& sublaw : mSublawIPs)
    {
        sublaw.ShiftToPast();
    }
}


void IPAdditiveInputImplicit::ShiftToFuture()
{
    for(auto& sublaw : mSublawIPs)
    {
        sublaw.ShiftToFuture();
    }
}


void IPAdditiveInputImplicit::NuToSerializeSave(SerializeStreamOut& rStream)
{
    IPConstitutiveLawBase::NuToSerializeSave(rStream);
    for(auto& sublaw : mSublawIPs)
    {
        sublaw.NuToSerializeSave(rStream);
    }
}


void IPAdditiveInputImplicit::NuToSerializeLoad(SerializeStreamIn& rStream)
{
    IPConstitutiveLawBase::NuToSerializeLoad(rStream);
    for(auto& sublaw : mSublawIPs)
    {
        sublaw.NuToSerializeLoad(rStream);
    }
}


template<int TDim>
void IPAdditiveInputImplicit::CalculateGlobalOutputs(const NuTo::ConstitutiveInputMap& rConstitutiveInput,
                                                     const NuTo::ConstitutiveOutputMap& rConstitutiveOutput,
                                                     std::vector<NuTo::ConstitutiveInputMap>& rLocalInputMapVec,
                                                     std::vector<ConstitutiveOutputMap>& rLocalOutputMapVec)
{
    static const bool SublawsHaveDamping = mLaw.CheckDofCombinationComputable(Node::eDof::DISPLACEMENTS,
                                                                              Node::eDof::DISPLACEMENTS,
                                                                              1);

    const constexpr int VoigtDim = ConstitutiveIOBase::GetVoigtDim(TDim);

    const int NumTimeDerivatives = SublawsHaveDamping ? 1:0;

    // Get time data
    double& currTime = mData.GetData(0).GetTime();
    const auto& itTime = rConstitutiveInput.find(Constitutive::eInput::TIME);
    if(itTime!=rConstitutiveInput.end())
        currTime = (*itTime->second)[0];
    double delta_t = currTime - mData.GetData(1).GetTime();



    // Get global outputs from local outputs
    for(const auto& itOutput : rConstitutiveOutput)
    {
        switch(itOutput.first)
        {
        case Constitutive::eOutput::ENGINEERING_STRESS:
        {
            const unsigned int numLocalUnknownsPerTimeDer = VoigtDim*mSublawIPs.size();
            unsigned int numLocalUnknowns = numLocalUnknownsPerTimeDer;

            if(rConstitutiveInput.find(eInput::CALCULATE_INITIALIZE_VALUE_RATES)!=rConstitutiveInput.end())
            {
                const unsigned int numTotalUnknowns = numLocalUnknowns + VoigtDim;

                // If static data vectors for local inputs have the wrong size the only reason should be, that they are not initialized!!
                // The problem is, that I have not found a way so far, to initialize them correctly before the first evaluate call.
                Eigen::VectorXd& SDCurrentStrain = mData.GetData(0).GetLocalInputs();
                Eigen::VectorXd& SDPreviousStrain = mData.GetData(1).GetLocalInputs();
                Eigen::VectorXd& SDCurrentStrainRate = mData.GetData(0).GetLocalInputRates();
                Eigen::VectorXd& SDPreviousStrainRate = mData.GetData(1).GetLocalInputRates();
                if(SDCurrentStrain.size()!=numLocalUnknownsPerTimeDer &&
                        SDPreviousStrain.size()!=numLocalUnknownsPerTimeDer)   // So far I have no idea, how to resize the vector before the first call of evaluate --- data during construction is missing
                {
                    SDCurrentStrain.resize(numLocalUnknownsPerTimeDer);
                    SDCurrentStrain.setZero();
                    SDPreviousStrain.resize(numLocalUnknownsPerTimeDer);
                    SDPreviousStrain.setZero();
                    SDCurrentStrainRate.resize(numLocalUnknownsPerTimeDer);
                    SDCurrentStrainRate.setZero();
                    SDPreviousStrainRate.resize(numLocalUnknownsPerTimeDer);
                    SDPreviousStrainRate.setZero();
                }


                // Generate lhs matrix
                // -------------------
                Eigen::MatrixXd lhsMat = Eigen::ArrayXXd::Zero(numTotalUnknowns, numTotalUnknowns);

                for (unsigned int i = 0; i < mSublawIPs.size(); ++i)
                {
                    unsigned int StartIndex = VoigtDim*i;


#ifndef NDEBUG
                    rLocalOutputMapVec[i].find(Constitutive::eOutput::D_ENGINEERING_STRESS_D_ENGINEERING_STRAIN_DT1)->second.get()->AssertIsMatrix<VoigtDim,VoigtDim>(Constitutive::eOutput::D_ENGINEERING_STRESS_D_ENGINEERING_STRAIN_DT1,
                                                                                                                                                                      __PRETTY_FUNCTION__);
#endif

                    lhsMat.block(StartIndex,StartIndex,VoigtDim,VoigtDim) = (*static_cast<ConstitutiveMatrix<VoigtDim,VoigtDim>*>(rLocalOutputMapVec[i].find(Constitutive::eOutput::D_ENGINEERING_STRESS_D_ENGINEERING_STRAIN_DT1)->second.get()));
                    lhsMat.block(numLocalUnknowns,StartIndex,VoigtDim,VoigtDim) = Eigen::MatrixXd::Identity(VoigtDim,VoigtDim);
                    lhsMat.block(StartIndex,numLocalUnknowns,VoigtDim,VoigtDim) = Eigen::MatrixXd::Identity(VoigtDim,VoigtDim) * -1.0;
                }




                // Generate rhs vector
                // -------------------
#ifndef NDEBUG
                // WATCHOUT There is no funtion which takes input enums, so I took the equivalent output enum even though the variable is an input!!!
                rConstitutiveInput.find(Constitutive::eInput::ENGINEERING_STRAIN_DT1)->second.get()->AssertIsVector<VoigtDim>(Constitutive::eOutput::ENGINEERING_STRAIN,
                                                                                                                              __PRETTY_FUNCTION__);
#endif
                const ConstitutiveVector<VoigtDim>& globalStrainRate = *static_cast<ConstitutiveVector<VoigtDim>*>(rConstitutiveInput.find(Constitutive::eInput::ENGINEERING_STRAIN_DT1)->second.get());
                Eigen::VectorXd rhsVec(numTotalUnknowns);
                rhsVec.setZero();


                for(unsigned int i=0; i<VoigtDim; ++i)
                {
                    rhsVec[numLocalUnknowns+i] = globalStrainRate[i];
                }

                for (unsigned int i = 0; i < mSublawIPs.size(); ++i)
                {
#ifndef NDEBUG
                    rLocalOutputMapVec[i].find(Constitutive::eOutput::D_ENGINEERING_STRESS_D_ENGINEERING_STRAIN)->second.get()->AssertIsMatrix<VoigtDim,VoigtDim>(Constitutive::eOutput::D_ENGINEERING_STRESS_D_ENGINEERING_STRAIN,
                                                                                                                                                                  __PRETTY_FUNCTION__);
#endif

                    Eigen::VectorXd localStrains(VoigtDim);
                    for(unsigned int j=0; j<VoigtDim; ++j)
                    {
                        localStrains[j] = SDCurrentStrain[i*VoigtDim +j];
                    }

                    auto localRHS = (*static_cast<ConstitutiveMatrix<VoigtDim,VoigtDim>*>(rLocalOutputMapVec[i].find(Constitutive::eOutput::D_ENGINEERING_STRESS_D_ENGINEERING_STRAIN_DT1)->second.get())) *
                            localStrains;

                    for(unsigned int j=0; j<VoigtDim; ++j)
                    {
                        rhsVec[i*VoigtDim + j] = localRHS[j];
                    }
                }

                for(unsigned int i=0; i<VoigtDim; ++i)
                {
                    rhsVec[numLocalUnknowns+i] = globalStrainRate[i];
                }



                // Solve local system
                // ------------------
                Eigen::PartialPivLU<Eigen::MatrixXd> Solver(lhsMat);
                Eigen::VectorXd resultVec =  Solver.solve(rhsVec);

                //            std::cout << lhsMat << std::endl <<"------------"<< std::endl
                //                      << rhsVec << std::endl <<"------------"<< std::endl
                //                      << resultVec << std::endl << std::endl << std::endl;

                Eigen::VectorXd& SDCurrentStresses = mData.GetData().GetStress();
                if(SDCurrentStresses.size()!=VoigtDim)
                    SDCurrentStresses.resize(VoigtDim);

                for (unsigned int i = 0; i < numLocalUnknownsPerTimeDer; ++i)
                {
                    mData.GetData(0).GetLocalInputRates()[i] = resultVec[i];
                    mData.GetData(1).GetLocalInputRates()[i] = resultVec[i];
                }

                for(unsigned int i=0; i<VoigtDim; ++i)
                {
                    (*itOutput.second)[i] = resultVec[i+numLocalUnknowns];
                    SDCurrentStresses[i] = resultVec[i+numLocalUnknowns];
                }
            }
            else
            {
                if (NumTimeDerivatives == 1 && delta_t >0.0)
                    numLocalUnknowns+=numLocalUnknownsPerTimeDer;
                const unsigned int numTotalUnknowns = numLocalUnknowns + VoigtDim;


                // Generate lhs matrix
                // -------------------
                Eigen::MatrixXd lhsMat = Eigen::ArrayXXd::Zero(numTotalUnknowns, numTotalUnknowns);



                for (unsigned int i = 0; i < mSublawIPs.size(); ++i)
                {
                    unsigned int StartIndex = VoigtDim*i;

#ifndef NDEBUG
                    rLocalOutputMapVec[i].find(Constitutive::eOutput::D_ENGINEERING_STRESS_D_ENGINEERING_STRAIN)->second.get()->AssertIsMatrix<VoigtDim,VoigtDim>(Constitutive::eOutput::D_ENGINEERING_STRESS_D_ENGINEERING_STRAIN,
                                                                                                                                                                  __PRETTY_FUNCTION__);
#endif
                    lhsMat.block(StartIndex,StartIndex,VoigtDim,VoigtDim) = (*static_cast<ConstitutiveMatrix<VoigtDim,VoigtDim>*>(rLocalOutputMapVec[i].find(Constitutive::eOutput::D_ENGINEERING_STRESS_D_ENGINEERING_STRAIN)->second.get()));
                    lhsMat.block(numLocalUnknowns,StartIndex,VoigtDim,VoigtDim) = Eigen::MatrixXd::Identity(VoigtDim,VoigtDim);
                    lhsMat.block(StartIndex,numLocalUnknowns,VoigtDim,VoigtDim) = Eigen::MatrixXd::Identity(VoigtDim,VoigtDim) * -1.0;
                    if(SublawsHaveDamping && delta_t >0)
                    {
                        unsigned int StartIndexFirstTimeDer = StartIndex+numLocalUnknownsPerTimeDer;
                        lhsMat.block(StartIndex,StartIndexFirstTimeDer,VoigtDim,VoigtDim) = (*static_cast<ConstitutiveMatrix<VoigtDim,VoigtDim>*>(rLocalOutputMapVec[i].find(Constitutive::eOutput::D_ENGINEERING_STRESS_D_ENGINEERING_STRAIN_DT1)->second.get()));
                    }
                }
                if(SublawsHaveDamping && delta_t>0)
                {
                    lhsMat.block(numLocalUnknownsPerTimeDer,0,numLocalUnknownsPerTimeDer,numLocalUnknownsPerTimeDer) = Eigen::MatrixXd::Identity(numLocalUnknownsPerTimeDer,numLocalUnknownsPerTimeDer);
                    lhsMat.block(numLocalUnknownsPerTimeDer,numLocalUnknownsPerTimeDer,numLocalUnknownsPerTimeDer,numLocalUnknownsPerTimeDer) = Eigen::MatrixXd::Identity(numLocalUnknownsPerTimeDer,numLocalUnknownsPerTimeDer) * (-delta_t);
                }


                // Generate rhs vector
                // -------------------
#ifndef NDEBUG
                // WATCHOUT There is no funtion which takes input enums, so I took the equivalent output enum even though the variable is an input!!!
                rConstitutiveInput.find(Constitutive::eInput::ENGINEERING_STRAIN)->second.get()->AssertIsVector<VoigtDim>(Constitutive::eOutput::ENGINEERING_STRAIN,
                                                                                                                          __PRETTY_FUNCTION__);
#endif
                const ConstitutiveVector<VoigtDim>& globalStrain = *static_cast<ConstitutiveVector<VoigtDim>*>(rConstitutiveInput.find(Constitutive::eInput::ENGINEERING_STRAIN)->second.get());
                Eigen::VectorXd rhsVec(numTotalUnknowns);
                rhsVec.setZero();

                for(unsigned int i=0; i<VoigtDim; ++i)
                {
                    rhsVec[numLocalUnknowns+i] = globalStrain[i];
                }
                if(SublawsHaveDamping)
                {
                    // If static data vectors for local inputs have the wrong size the only reason should be, that they are not initialized!!
                    // The problem is, that I have not found a way so far, to initialize them correctly before the first evaluate call.
                    Eigen::VectorXd& SDCurrentStrain = mData.GetData(0).GetLocalInputs();
                    Eigen::VectorXd& SDPreviousStrain = mData.GetData(1).GetLocalInputs();
                    Eigen::VectorXd& SDCurrentStrainRate = mData.GetData(0).GetLocalInputRates();
                    Eigen::VectorXd& SDPreviousStrainRate = mData.GetData(1).GetLocalInputRates();
                    if(SDCurrentStrain.size()!=numLocalUnknownsPerTimeDer &&
                            SDPreviousStrain.size()!=numLocalUnknownsPerTimeDer)   // So far I have no idea, how to resize the vector before the first call of evaluate --- data during construction is missing
                    {
                        SDCurrentStrain.resize(numLocalUnknownsPerTimeDer);
                        SDCurrentStrain.setZero();
                        SDPreviousStrain.resize(numLocalUnknownsPerTimeDer);
                        SDPreviousStrain.setZero();
                        SDCurrentStrainRate.resize(numLocalUnknownsPerTimeDer);
                        SDCurrentStrainRate.setZero();
                        SDPreviousStrainRate.resize(numLocalUnknownsPerTimeDer);
                        SDPreviousStrainRate.setZero();
                    }
                    if(delta_t > 0)
                    {
                        for (unsigned int i = 0; i < numLocalUnknownsPerTimeDer; ++i)
                        {
                            rhsVec[i+numLocalUnknownsPerTimeDer] = SDPreviousStrain[i];
                        }
                    }
                    else
                    {
                        for (unsigned int i = 0; i < mSublawIPs.size(); ++i)
                        {
#ifndef NDEBUG
                            rLocalOutputMapVec[i].find(Constitutive::eOutput::D_ENGINEERING_STRESS_D_ENGINEERING_STRAIN_DT1)->second.get()->AssertIsMatrix<VoigtDim,VoigtDim>(Constitutive::eOutput::D_ENGINEERING_STRESS_D_ENGINEERING_STRAIN_DT1,
                                                                                                                                                                          __PRETTY_FUNCTION__);
#endif

                            SDCurrentStrainRate = SDPreviousStrainRate;

                            Eigen::VectorXd localStrainRates(VoigtDim);
                            for(unsigned int j=0; j<VoigtDim; ++j)
                            {
                                localStrainRates[j] = SDCurrentStrainRate[i*VoigtDim +j];
                            }

                            Eigen::VectorXd localRHS = (*static_cast<ConstitutiveMatrix<VoigtDim,VoigtDim>*>(rLocalOutputMapVec[i].find(Constitutive::eOutput::D_ENGINEERING_STRESS_D_ENGINEERING_STRAIN_DT1)->second.get())) *
                                             localStrainRates *-1.0;

                            for(unsigned int j=0; j<VoigtDim; ++j)
                            {
                                rhsVec[i*VoigtDim + j] = localRHS[j];
                            }
                        }
                    }
                }

                // Solve local system
                // ------------------

                //            Eigen::VectorXd resultVec = lhsMat.inverse() * rhsVec;

                Eigen::PartialPivLU<Eigen::MatrixXd> Solver(lhsMat);
                Eigen::VectorXd resultVec =  Solver.solve(rhsVec);

//                            std::cout << lhsMat << std::endl <<"------------"<< std::endl
//                                      << rhsVec << std::endl <<"------------"<< std::endl
//                                      << resultVec << std::endl << std::endl << std::endl;

                Eigen::VectorXd& SDCurrentStresses = mData.GetData().GetStress();
                if(SDCurrentStresses.size()!=VoigtDim)
                    SDCurrentStresses.resize(VoigtDim);
                if(SublawsHaveDamping && delta_t > 0)
                {
                    for (unsigned int i = 0; i < numLocalUnknownsPerTimeDer; ++i)
                    {
                        mData.GetData(0).GetLocalInputs()[i] = resultVec[i];
                        mData.GetData(0).GetLocalInputRates()[i] = resultVec[i +numLocalUnknownsPerTimeDer];
                    }
                }
                for(unsigned int i=0; i<VoigtDim; ++i)
                {
                    (*itOutput.second)[i] = resultVec[i+numLocalUnknowns];
                    SDCurrentStresses[i] = resultVec[i+numLocalUnknowns];
                }
            }
            itOutput.second->SetIsCalculated(true);
            break;
        }


        case Constitutive::eOutput::ENGINEERING_STRESS_VISUALIZE:
        {
            Eigen::VectorXd& SDCurrentStresses = mData.GetData().GetStress();

            switch(TDim)
            {
            case 1:
                (*itOutput.second)[0] = SDCurrentStresses[0];
                (*itOutput.second)[1] = 0;
                (*itOutput.second)[2] = 0;
                (*itOutput.second)[3] = 0;
                (*itOutput.second)[4] = 0;
                (*itOutput.second)[5] = 0;
                break;

            case 2:
                (*itOutput.second)[0] = SDCurrentStresses[0];
                (*itOutput.second)[1] = SDCurrentStresses[1];
                (*itOutput.second)[2] = 0;
                (*itOutput.second)[3] = SDCurrentStresses[2];
                (*itOutput.second)[4] = 0;
                (*itOutput.second)[5] = 0;
                break;

            case 3:
                for(unsigned int i=0; i<VoigtDim; ++i)
                {
                    (*itOutput.second)[i] = SDCurrentStresses[i];
                }
                break;

            default:
                break;
            }
            itOutput.second->SetIsCalculated(true);
            break;
        }


        case Constitutive::eOutput::D_ENGINEERING_STRESS_D_ENGINEERING_STRAIN:
        {

            Eigen::Matrix<double,VoigtDim,VoigtDim> globalCompliance(VoigtDim,VoigtDim);
            globalCompliance.setZero();
            for(const auto& localOutputs : rLocalOutputMapVec)
            {
                Eigen::Matrix<double,VoigtDim,VoigtDim>& localStiffness = (*static_cast<ConstitutiveMatrix<VoigtDim,VoigtDim>*>(localOutputs.find(Constitutive::eOutput::D_ENGINEERING_STRESS_D_ENGINEERING_STRAIN)->second.get())).AsMatrix();

                if(SublawsHaveDamping && delta_t >0)
                {
                    Eigen::Matrix<double,VoigtDim,VoigtDim>& localDamping = (*static_cast<ConstitutiveMatrix<VoigtDim,VoigtDim>*>(localOutputs.find(Constitutive::eOutput::D_ENGINEERING_STRESS_D_ENGINEERING_STRAIN_DT1)->second.get())).AsMatrix();
                    localStiffness += localDamping / delta_t;
                }
                globalCompliance += localStiffness.inverse();
            }

            itOutput.second->AssertIsMatrix<VoigtDim,VoigtDim>(itOutput.first,__PRETTY_FUNCTION__);
            Eigen::Matrix<double,VoigtDim,VoigtDim>& outputMat = (*static_cast<ConstitutiveMatrix<VoigtDim,VoigtDim>*>(itOutput.second.get())).AsMatrix();

            outputMat = globalCompliance.inverse();
            itOutput.second->SetIsCalculated(true);
            break;
        }

        case Constitutive::eOutput::D_ENGINEERING_STRESS_D_ENGINEERING_STRAIN_DT1:
        {
            itOutput.second->AssertIsMatrix<VoigtDim,VoigtDim>(itOutput.first,__PRETTY_FUNCTION__);
            Eigen::Matrix<double,VoigtDim,VoigtDim>& outputMat = (*static_cast<ConstitutiveMatrix<VoigtDim,VoigtDim>*>(itOutput.second.get())).AsMatrix();
            if(rConstitutiveInput.find(eInput::CALCULATE_INITIALIZE_VALUE_RATES)!=rConstitutiveInput.end())
            {
                Eigen::Matrix<double,VoigtDim,VoigtDim> globalCompliance(VoigtDim,VoigtDim);
                globalCompliance.setZero();
                for(const auto& localOutputs : rLocalOutputMapVec)
                {

                    Eigen::Matrix<double,VoigtDim,VoigtDim>& localDamping = (*static_cast<ConstitutiveMatrix<VoigtDim,VoigtDim>*>(localOutputs.find(Constitutive::eOutput::D_ENGINEERING_STRESS_D_ENGINEERING_STRAIN_DT1)->second.get())).AsMatrix();

                    globalCompliance += localDamping.inverse();
                }
                outputMat = globalCompliance.inverse();
                itOutput.second->SetIsCalculated(true);
                break;
            }
            outputMat.setZero();
            itOutput.second->SetIsCalculated(true);
            break;
        }

        case Constitutive::eOutput::UPDATE_STATIC_DATA:
        {
            mData.GetData(1) = mData.GetData();
        }

        default:
            break;
        }
    }


}


template<int TDim>
void IPAdditiveInputImplicit::CreateLocalInAndOutputMaps(const NuTo::ConstitutiveInputMap& rConstitutiveInput,
                                                         const NuTo::ConstitutiveOutputMap& rConstitutiveOutput,
                                                         std::vector<NuTo::ConstitutiveInputMap>& rLocalInputMapVec,
                                                         std::vector<ConstitutiveOutputMap>& rLocalOutputMapVec)
{
    static const bool SublawsHaveDamping = mLaw.CheckDofCombinationComputable(Node::eDof::DISPLACEMENTS,
                                                                              Node::eDof::DISPLACEMENTS,
                                                                              1);

    const constexpr int VoigtDim = ConstitutiveIOBase::GetVoigtDim(TDim);
    assert(rLocalInputMapVec.size() ==mSublawIPs.size() && "Referenced vector should already have the necessary size (Resize during construction!)");
    assert(rLocalOutputMapVec.size()==mSublawIPs.size() && "Referenced vector should already have the necessary size (Resize during construction!)");

    // Copy inputs for every constitutive law
    for (unsigned int i = 0; i < mSublawIPs.size(); ++i)
    {
        for (const auto& itInput : rConstitutiveInput)
        {
            if(itInput.first == Constitutive::eInput::CALCULATE_INITIALIZE_VALUE_RATES)
                rLocalInputMapVec[i].emplace(itInput.first, nullptr);
            else
                rLocalInputMapVec[i].emplace(itInput.first, itInput.second->clone());
        }
    }

    // Copy outputs for every constitutive law and add additional necessary outputs
    for (unsigned int i = 0; i < mSublawIPs.size(); ++i)
    {
        for(const auto& itOutput: rConstitutiveOutput)
        {
            switch(itOutput.first)
            {
            case Constitutive::eOutput::UPDATE_STATIC_DATA:
                rLocalOutputMapVec[i].emplace(itOutput.first, nullptr);
                break;

            case Constitutive::eOutput::ENGINEERING_STRESS:
                rLocalOutputMapVec[i].emplace(itOutput.first, itOutput.second->clone());
                rLocalOutputMapVec[i].emplace(Constitutive::eOutput::D_ENGINEERING_STRESS_D_ENGINEERING_STRAIN,
                                              std::make_unique<ConstitutiveMatrix<VoigtDim,VoigtDim>>());
                if(SublawsHaveDamping)
                    rLocalOutputMapVec[i].emplace(Constitutive::eOutput::D_ENGINEERING_STRESS_D_ENGINEERING_STRAIN_DT1,
                                                  std::make_unique<ConstitutiveMatrix<VoigtDim,VoigtDim>>());

                break;


                // no need for local copy!
            case Constitutive::eOutput::ENGINEERING_STRESS_VISUALIZE:
                break;

            default:
                rLocalOutputMapVec[i].emplace(itOutput.first, itOutput.second->clone());
                break;
            }
        }
    }
}




template<int TDim>
NuTo::eError IPAdditiveInputImplicit::AdditiveInputImplicitEvaluate(
        const NuTo::ConstitutiveInputMap& rConstitutiveInput,
        const NuTo::ConstitutiveOutputMap& rConstitutiveOutput)
{
    static_assert (TDim == 1 || TDim == 2 || TDim == 3 , "Dimensions 1D, 2D & 3D supported.");
    NuTo::eError error = NuTo::eError::SUCCESSFUL;
//    const constexpr int VoigtDim = ConstitutiveIOBase::GetVoigtDim(TDim);


    std::vector<ConstitutiveInputMap>  localInputMapVec(mSublawIPs.size());
    std::vector<ConstitutiveOutputMap> localOutputMapVec(mSublawIPs.size());

    // Get local inputs and outputs
    CreateLocalInAndOutputMaps<TDim>(rConstitutiveInput,
                                     rConstitutiveOutput,
                                     localInputMapVec,
                                     localOutputMapVec);




    // evaluate sublaws
    for (unsigned int i = 0; i < mSublawIPs.size(); ++i)
    {
        eError error = mSublawIPs[i].Evaluate<TDim>(localInputMapVec[i], localOutputMapVec[i]);
        if(error!=eError::SUCCESSFUL)
            throw Exception(__PRETTY_FUNCTION__,
                            "One or more attached constitutive laws return error codes. Can't handle this");
    }


    // calculate global outputs
    CalculateGlobalOutputs<TDim>(rConstitutiveInput,
                                 rConstitutiveOutput,
                                 localInputMapVec,
                                 localOutputMapVec);



    return error;
}







template NuTo::eError IPAdditiveInputImplicit::AdditiveInputImplicitEvaluate<1>(
        const NuTo::ConstitutiveInputMap& rConstitutiveInput,
        const NuTo::ConstitutiveOutputMap& rConstitutiveOutput);
template NuTo::eError IPAdditiveInputImplicit::AdditiveInputImplicitEvaluate<2>(
        const NuTo::ConstitutiveInputMap& rConstitutiveInput,
        const NuTo::ConstitutiveOutputMap& rConstitutiveOutput);
template NuTo::eError IPAdditiveInputImplicit::AdditiveInputImplicitEvaluate<3>(
        const NuTo::ConstitutiveInputMap& rConstitutiveInput,
        const NuTo::ConstitutiveOutputMap& rConstitutiveOutput);
