#include "AdditiveInputImplicit.h"
#include "nuto/mechanics/constitutive/inputoutput/EngineeringStrain.h"
#include "nuto/mechanics/constitutive/inputoutput/EngineeringStress.h"


#include "nuto/math/SparseDirectSolverMUMPS.h"
#include "nuto/math/SparseDirectSolverMKLPardiso.h"
#include "nuto/math/SparseDirectSolverPardiso.h"
#include "nuto/math/SparseMatrixCSRGeneral.h"


void NuTo::AdditiveInputImplicit::AddConstitutiveLaw(NuTo::ConstitutiveBase *rConstitutiveLaw, Constitutive::Input::eInput rModiesInput)
{
    if(mStaticDataAllocated)
        throw MechanicsException(__PRETTY_FUNCTION__,"All constitutive laws have to be attached before static data is allocated!");
    mConstitutiveLaws.push_back(rConstitutiveLaw);
    AddCalculableDofCombinations(rConstitutiveLaw);
}


bool NuTo::AdditiveInputImplicit::CheckDofCombinationComputable(NuTo::Node::eDof rDofRow, NuTo::Node::eDof rDofCol, int rTimeDerivative) const
{
    if(mComputableDofCombinations[rTimeDerivative].find(std::pair<Node::eDof,Node::eDof>(rDofRow,rDofCol)) != mComputableDofCombinations[rTimeDerivative].end())
        return true;
    return false;
}

template <int TDim>
NuTo::Error::eError NuTo::AdditiveInputImplicit::EvaluateAdditiveInputImplicit(NuTo::ElementBase *rElement,
                                                                               int rIp,
                                                                               const NuTo::ConstitutiveInputMap &rConstitutiveInput,
                                                                               const NuTo::ConstitutiveOutputMap &rConstitutiveOutput)
{
    static_assert (TDim == 1 || TDim == 2 || TDim == 3 , "Dimensions 1D, 2D & 3D supported.");

    Error::eError error = Error::SUCCESSFUL;
    const constexpr int VoigtDim = ConstitutiveIOBase::GetVoigtDim(TDim);


    // Copy inputs for every constitutive law
    std::vector<NuTo::ConstitutiveInputMap> localInputMapVec(mConstitutiveLaws.size());
    for(unsigned int i=0; i<mConstitutiveLaws.size(); ++i)
    {
        for(const auto& itInput: rConstitutiveInput)
        {
            localInputMapVec[i].emplace(itInput.first,
                                       (*itInput.second).clone());
        }
    }



    // Copy outputs for every constitutive law
    std::vector<ConstitutiveOutputMap> localOutputMapVec(mConstitutiveLaws.size());
    for(unsigned int i=0; i<mConstitutiveLaws.size(); ++i)
    {
        for(const auto& itOutput: rConstitutiveOutput)
        {
            if(itOutput.first == Constitutive::Output::UPDATE_STATIC_DATA)
            {
                localOutputMapVec[i].emplace(itOutput.first,
                                             nullptr);
            }
            else
            {
                localOutputMapVec[i].emplace(itOutput.first,
                                             (*itOutput.second).clone());
            }

            // temporary: until a better solution is found
            if(itOutput.first == Constitutive::Output::ENGINEERING_STRESS &&
               rConstitutiveOutput.find(Constitutive::Output::D_ENGINEERING_STRESS_D_ENGINEERING_STRAIN) == rConstitutiveOutput.end())
                localOutputMapVec[i].emplace(Constitutive::Output::D_ENGINEERING_STRESS_D_ENGINEERING_STRAIN,
                                             std::make_unique<ConstitutiveMatrix<VoigtDim,VoigtDim>>());

            // temporary: Implementation of static data not completed
            if(itOutput.first == Constitutive::Output::D_ENGINEERING_STRESS_D_ENGINEERING_STRAIN &&
               rConstitutiveOutput.find(Constitutive::Output::ENGINEERING_STRESS) == rConstitutiveOutput.end())
                localOutputMapVec[i].emplace(Constitutive::Output::ENGINEERING_STRESS,
                                             std::make_unique<ConstitutiveVector<VoigtDim>>());


            // temporary: Implementation of static data not completed
            if(itOutput.first == Constitutive::Output::ENGINEERING_STRESS_VISUALIZE)
            {
                if(rConstitutiveOutput.find(Constitutive::Output::ENGINEERING_STRESS) == rConstitutiveOutput.end())
                    localOutputMapVec[i].emplace(Constitutive::Output::ENGINEERING_STRESS,
                                                 std::make_unique<ConstitutiveVector<VoigtDim>>());


                if (rConstitutiveOutput.find(Constitutive::Output::D_ENGINEERING_STRESS_D_ENGINEERING_STRAIN) == rConstitutiveOutput.end())
                    localOutputMapVec[i].emplace(Constitutive::Output::D_ENGINEERING_STRESS_D_ENGINEERING_STRAIN,
                                                 std::make_unique<ConstitutiveMatrix<VoigtDim,VoigtDim>>());
            }
        }


    }


    for(unsigned int i=0; i<mConstitutiveLaws.size(); ++i)
    {
        error = mConstitutiveLaws[i]->Evaluate<TDim>(rElement,
                                                     rIp,
                                                     localInputMapVec[i],
                                                     localOutputMapVec[i]);
        if(error!=Error::SUCCESSFUL)
            throw Exception(__PRETTY_FUNCTION__,
                            "One or more attached constitutive laws return error codes. Can't handle this");
    }

    if(rConstitutiveOutput.find(Constitutive::Output::ENGINEERING_STRESS) != rConstitutiveOutput.end() ||
       rConstitutiveOutput.find(Constitutive::Output::ENGINEERING_STRESS_VISUALIZE) != rConstitutiveOutput.end() ||
       rConstitutiveOutput.find(Constitutive::Output::D_ENGINEERING_STRESS_D_ENGINEERING_STRAIN) != rConstitutiveOutput.end())
    {

        double numLocalUnknowns = VoigtDim*mConstitutiveLaws.size();
        double numTotalUnknowns = numLocalUnknowns + VoigtDim;


        // Generate lhs matrix
        Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> lhsMat = Eigen::ArrayXXd::Zero(numTotalUnknowns, numTotalUnknowns);

        for(unsigned int i= 0; i<mConstitutiveLaws.size(); ++i)
        {
            unsigned int StartIndex = VoigtDim*i;

#ifdef DEBUG
            localOutputMapVec[i].find(Constitutive::Output::D_ENGINEERING_STRESS_D_ENGINEERING_STRAIN)->second.get()->AssertIsMatrix<VoigtDim,VoigtDim>(Constitutive::Output::D_ENGINEERING_STRESS_D_ENGINEERING_STRAIN,
                                                                                                                                                        __PRETTY_FUNCTION__);
#endif
            lhsMat.block(StartIndex,StartIndex,VoigtDim,VoigtDim) = (*static_cast<ConstitutiveMatrix<VoigtDim,VoigtDim>*>(localOutputMapVec[i].find(Constitutive::Output::D_ENGINEERING_STRESS_D_ENGINEERING_STRAIN)->second.get()));
            lhsMat.block(numLocalUnknowns,StartIndex,VoigtDim,VoigtDim) = Eigen::MatrixXd::Identity(VoigtDim,VoigtDim);
            lhsMat.block(StartIndex,numLocalUnknowns,VoigtDim,VoigtDim) = Eigen::MatrixXd::Identity(VoigtDim,VoigtDim) * -1.0;
        }


        // Generate rhs vector
#ifdef DEBUG
        // WATCHOUT There is no funtion which takes input enums, so I took the equivalent output enum even though the variable is an input!!!
        rConstitutiveInput.find(Constitutive::Input::ENGINEERING_STRAIN)->second.get()->AssertIsVector<VoigtDim>(Constitutive::Output::ENGINEERING_STRAIN,
                                                                                                                 __PRETTY_FUNCTION__);
#endif
        const ConstitutiveVector<VoigtDim>& globalStrain = *static_cast<ConstitutiveVector<VoigtDim>*>(rConstitutiveInput.find(Constitutive::Input::ENGINEERING_STRAIN)->second.get());
        FullVector<double,Eigen::Dynamic> rhsVec(numTotalUnknowns);
        rhsVec.setZero();
        for(unsigned int i=0; i<VoigtDim; ++i)
        {
            rhsVec[numLocalUnknowns+i] = globalStrain[i];
        }

        FullVector<double,Eigen::Dynamic> resultVec;


#if defined(HAVE_PARDISO) && defined(_OPENMP)
        const unsigned int NumProcessors = 4;
        NuTo::SparseDirectSolverPardiso mySolver(NumProcessors, 0); // note: not the MKL version
        mySolver.SetShowTime(false);
#else
        NuTo::SparseDirectSolverMUMPS mySolver;
        mySolver.SetVerboseLevel(0);
        mySolver.SetShowTime(false);
#endif

        NuTo::SparseMatrixCSRGeneral<double> matForSolver(lhsMat);
        matForSolver.SetOneBasedIndexing();

        //    auto test = lhsMat.inverse();
        //    FullVector<double,Eigen::Dynamic> result2 = test * rhsVec;

        mySolver.Solve(matForSolver,rhsVec,resultVec);



        for(const auto& itOutput : rConstitutiveOutput)
        {
            switch(itOutput.first)
            {
            case Constitutive::Output::ENGINEERING_STRESS:

                for(unsigned int i=0; i<VoigtDim; ++i)
                {
                    (*itOutput.second)[i] = resultVec[i+numLocalUnknowns];
                }
                itOutput.second->SetIsCalculated(true);
                break;

            case Constitutive::Output::D_ENGINEERING_STRESS_D_ENGINEERING_STRAIN:
            {

                Eigen::Matrix<double,VoigtDim,VoigtDim> BlockProduct(VoigtDim,VoigtDim);
                BlockProduct = -1 * lhsMat.block(numLocalUnknowns,0,VoigtDim,numLocalUnknowns)
                        * (lhsMat.block(0,0,numLocalUnknowns,numLocalUnknowns)).inverse()
                        * lhsMat.block(0,numLocalUnknowns,numLocalUnknowns,VoigtDim);

                itOutput.second->AssertIsMatrix<VoigtDim,VoigtDim>(itOutput.first,__PRETTY_FUNCTION__);
                Eigen::Matrix<double,VoigtDim,VoigtDim>& outputMat = (*static_cast<ConstitutiveMatrix<VoigtDim,VoigtDim>*>(itOutput.second.get())).AsMatrix();

                outputMat = BlockProduct.inverse();
                itOutput.second->SetIsCalculated(true);
            }
                break;

            case Constitutive::Output::ENGINEERING_STRESS_VISUALIZE:

                switch(TDim)
                {
                case 1:
                    (*itOutput.second)[0] = resultVec[numLocalUnknowns];
                    (*itOutput.second)[1] = 0;
                    (*itOutput.second)[2] = 0;
                    (*itOutput.second)[3] = 0;
                    (*itOutput.second)[4] = 0;
                    (*itOutput.second)[5] = 0;
                    break;

                case 2:
                    (*itOutput.second)[0] = resultVec[numLocalUnknowns];
                    (*itOutput.second)[1] = resultVec[numLocalUnknowns+1];
                    (*itOutput.second)[2] = 0;
                    (*itOutput.second)[3] = resultVec[numLocalUnknowns+2];
                    (*itOutput.second)[4] = 0;
                    (*itOutput.second)[5] = 0;
                    break;

                case 3:
                    for(unsigned int i=0; i<VoigtDim; ++i)
                    {
                        (*itOutput.second)[i] = resultVec[i+numLocalUnknowns];
                    }
                    break;

                default:
                    break;
                }
                itOutput.second->SetIsCalculated(true);
                break;


            default:
                break;
            }

        }
    }
    return error;
    // static data -> save previous values.

}




NuTo::ConstitutiveInputMap NuTo::AdditiveInputImplicit::GetConstitutiveInputs(const NuTo::ConstitutiveOutputMap &rConstitutiveOutput,
                                                                              const NuTo::InterpolationType &rInterpolationType) const
{
    ConstitutiveInputMap constitutiveInputMap;

    for(unsigned int i=0; i<mConstitutiveLaws.size(); ++i)
    {

        ConstitutiveInputMap singleLawInputMap = mConstitutiveLaws[i]->GetConstitutiveInputs(rConstitutiveOutput,
                                                                                             rInterpolationType);

        constitutiveInputMap.insert(std::move_iterator<ConstitutiveInputMap::iterator>(singleLawInputMap.begin()),
                                    std::move_iterator<ConstitutiveInputMap::iterator>(singleLawInputMap.end()));
    }
//    // VHIRTHAMTODO: Remove! --- Temporary because static data not implemented correct. Therefore iteration for local values must e performed during each evaluation!
    const auto& itDEngineeringStressDEngineeringStrain = rConstitutiveOutput.find(Constitutive::Output::D_ENGINEERING_STRESS_D_ENGINEERING_STRAIN);
    if(itDEngineeringStressDEngineeringStrain != rConstitutiveOutput.end())
    {
        constitutiveInputMap.emplace(Constitutive::Input::ENGINEERING_STRAIN,nullptr);
    }



    return constitutiveInputMap;
}




void NuTo::AdditiveInputImplicit::AddCalculableDofCombinations(NuTo::ConstitutiveBase *rConstitutiveLaw)
{
    std::set<Node::eDof> allDofs = Node::GetDofSet();
    for (unsigned int i=0; i<mComputableDofCombinations.size(); ++i)
    for (auto itRow : allDofs)
        for (auto itCol : allDofs)
        {
            if (rConstitutiveLaw->CheckDofCombinationComputable(itRow,itCol,i))
                    mComputableDofCombinations[i].emplace(itRow,itCol);
        }
}
