#pragma once


#include "nuto/base/Timer.h"
#include "nuto/mechanics/MechanicsException.h"
#include "nuto/mechanics/constitutive/inputoutput/ConstitutiveScalar.h"
#include "nuto/mechanics/constitutive/inputoutput/ConstitutiveVector.h"
#include "nuto/mechanics/constitutive/staticData/ConstitutiveStaticDataMoistureTransport.h"
#include "nuto/mechanics/elements/ElementOutputBlockVectorDouble.h"
#include "nuto/mechanics/elements/ElementOutputBlockMatrixDouble.h"
#include "nuto/mechanics/structures/unstructured/Structure.h"
#include "ConstitutiveLawMoistureTransport_Setup.h"










template<int TDim>
struct ConstitutiveOutputTest_EvaluateData
{
    //input
    NuTo::ConstitutiveScalar          mRelativeHumidity;
    NuTo::ConstitutiveScalar          mRelativeHumidity_dt1;
    NuTo::ConstitutiveVector<TDim>    mRelativeHumidity_Gradient;
    NuTo::ConstitutiveScalar          mWaterVolumeFraction;
    NuTo::ConstitutiveScalar          mWaterVolumeFraction_dt1;
    NuTo::ConstitutiveVector<TDim>    mWaterVolumeFraction_Gradient;
    //output - internal gradient
    NuTo::ConstitutiveVector<TDim>    mInternalGradientRH_B;
    NuTo::ConstitutiveScalar          mInternalGradientRH_N;
    NuTo::ConstitutiveVector<TDim>    mInternalGradientWV_B;
    NuTo::ConstitutiveScalar          mInternalGradientWV_N;
    //output - hessian 0
    NuTo::ConstitutiveScalar          mInternalGradientRH_dRH_BB_H0;
    NuTo::ConstitutiveScalar          mInternalGradientRH_dRH_NN_H0;
    NuTo::ConstitutiveVector<TDim>    mInternalGradientRH_dWV_BN_H0;
    NuTo::ConstitutiveScalar          mInternalGradientRH_dWV_NN_H0;
    NuTo::ConstitutiveScalar          mInternalGradientWV_dRH_NN_H0;
    NuTo::ConstitutiveScalar          mInternalGradientWV_dWV_BB_H0;
    NuTo::ConstitutiveVector<TDim>    mInternalGradientWV_dWV_BN_H0;
    NuTo::ConstitutiveScalar          mInternalGradientWV_dWV_NN_H0;
    //output - hessian 1
    NuTo::ConstitutiveScalar          mInternalGradientRH_dRH_NN_H1;
    NuTo::ConstitutiveScalar          mInternalGradientRH_dWV_NN_H1;
    NuTo::ConstitutiveScalar          mInternalGradientWV_dWV_NN_H1;
    //output - internal gradient - boundary
    NuTo::ConstitutiveScalar          mInternalGradientRH_Boundary_N;
    NuTo::ConstitutiveScalar          mInternalGradientWV_Boundary_N;
    //output - hessian 0 - boundary
    NuTo::ConstitutiveScalar          mInternalGradientRH_dRH_Boundary_NN_H0;
    NuTo::ConstitutiveScalar          mInternalGradientWV_dWV_Boundary_NN_H0;
};










void ConstitutiveOutputTest_SetupMoistureTransport(MoistureTransportControl& rMT)
{

    rMT.AdsorptionCoeffs = Eigen::MatrixXd::Constant(4, 1, 0.0);
    rMT.AdsorptionCoeffs(1) = 0.5;
    rMT.DesorptionCoeffs = Eigen::MatrixXd::Constant(4, 1, 0.0);
    rMT.DesorptionCoeffs(1) = 0.5;


    rMT.InitialRelativeHumidity = 1.0;
    rMT.MassExchangeRate = 1e-7;
    rMT.PoreVolumeFraction = 0.25;
    rMT.DiffusionCoefficientRH = 1e-12;
    rMT.DiffusionExponentRH = 1.0;
    rMT.DensitySaturatedWaterVapor = 0.2;
    rMT.DensityWater = 1000.0;
    rMT.DiffusionCoefficientWV = 1e-7;
    rMT.DiffusionExponentWV = 2.0;
    rMT.BoundaryEnvironmentalRH = 0.4;
    rMT.BoundaryDiffusionCoefficientRH = 1e-3;
    rMT.BoundaryDiffusionCoefficientWV = 1e-9;

    rMT.SetParametersConstitutiveLaw();
}


void ConstitutiveOutputTest_SetNodalValues(NuTo::Structure& rS, const MoistureTransportControl& rMT)
{
    for (int i=0; i<rS.GetNumNodes(); i++)
    {

        NuTo::NodeBase& node =  *rS.NodeGetNodePtr(i);

        double factor = (node.GetCoordinate(0)<0)? 0.5 : 1.0;

        if(node.GetNumRelativeHumidity() != 0)
        {
            node.SetRelativeHumidity(0,rMT.InitialRelativeHumidity * factor) ;
            node.SetRelativeHumidity(1,rMT.InitialRelativeHumidity * factor * 0.01) ;
        }
        if(node.GetNumWaterVolumeFraction() != 0)
        {
            node.SetWaterVolumeFraction(0,rMT.InitialWaterVolumeFraction * factor);
            node.SetWaterVolumeFraction(1,rMT.InitialWaterVolumeFraction * factor * 0.01);
        }
    }
}

template<int TDim>
NuTo::ConstitutiveInputMap ConstitutiveOutputTest_CreateConstitutiveInputMap(ConstitutiveOutputTest_EvaluateData<TDim>& evaluateData)
{
    NuTo::ConstitutiveInputMap constitutiveInputMap;
    constitutiveInputMap[NuTo::Constitutive::Input::RELATIVE_HUMIDITY]              = &(evaluateData.mRelativeHumidity);
    constitutiveInputMap[NuTo::Constitutive::Input::RELATIVE_HUMIDITY_DT1]          = &(evaluateData.mRelativeHumidity_dt1);
    constitutiveInputMap[NuTo::Constitutive::Input::RELATIVE_HUMIDITY_GRADIENT]     = &(evaluateData.mRelativeHumidity_Gradient);
    constitutiveInputMap[NuTo::Constitutive::Input::WATER_VOLUME_FRACTION]          = &(evaluateData.mWaterVolumeFraction);
    constitutiveInputMap[NuTo::Constitutive::Input::WATER_VOLUME_FRACTION_DT1]      = &(evaluateData.mWaterVolumeFraction_dt1);
    constitutiveInputMap[NuTo::Constitutive::Input::WATER_VOLUME_FRACTION_GRADIENT] = &(evaluateData.mWaterVolumeFraction_Gradient);
    return constitutiveInputMap;
}

template<int TDim>
NuTo::ConstitutiveOutputMap ConstitutiveOutputTest_CreateConstitutiveOutputMap(ConstitutiveOutputTest_EvaluateData<TDim>& evaluateData)
{
    NuTo::ConstitutiveOutputMap constitutiveOutputMap;
    constitutiveOutputMap[NuTo::Constitutive::Output::INTERNAL_GRADIENT_RELATIVE_HUMIDITY_B]              = &evaluateData.mInternalGradientRH_B;
    constitutiveOutputMap[NuTo::Constitutive::Output::INTERNAL_GRADIENT_RELATIVE_HUMIDITY_N]              = &evaluateData.mInternalGradientRH_N;
    constitutiveOutputMap[NuTo::Constitutive::Output::INTERNAL_GRADIENT_WATER_VOLUME_FRACTION_B]          = &evaluateData.mInternalGradientWV_B;
    constitutiveOutputMap[NuTo::Constitutive::Output::INTERNAL_GRADIENT_WATER_VOLUME_FRACTION_N]          = &evaluateData.mInternalGradientWV_N;

    constitutiveOutputMap[NuTo::Constitutive::Output::D_INTERNAL_GRADIENT_RH_D_RH_BB_H0]                  = &evaluateData.mInternalGradientRH_dRH_BB_H0;
    constitutiveOutputMap[NuTo::Constitutive::Output::D_INTERNAL_GRADIENT_RH_D_RH_NN_H0]                  = &evaluateData.mInternalGradientRH_dRH_NN_H0;
    constitutiveOutputMap[NuTo::Constitutive::Output::D_INTERNAL_GRADIENT_RH_D_WV_BN_H0]                  = &evaluateData.mInternalGradientRH_dWV_BN_H0;
    constitutiveOutputMap[NuTo::Constitutive::Output::D_INTERNAL_GRADIENT_RH_D_WV_NN_H0]                  = &evaluateData.mInternalGradientRH_dWV_NN_H0;
    constitutiveOutputMap[NuTo::Constitutive::Output::D_INTERNAL_GRADIENT_WV_D_RH_NN_H0]                  = &evaluateData.mInternalGradientWV_dRH_NN_H0;
    constitutiveOutputMap[NuTo::Constitutive::Output::D_INTERNAL_GRADIENT_WV_D_WV_BB_H0]                  = &evaluateData.mInternalGradientWV_dWV_BB_H0;
    constitutiveOutputMap[NuTo::Constitutive::Output::D_INTERNAL_GRADIENT_WV_D_WV_BN_H0]                  = &evaluateData.mInternalGradientWV_dWV_BN_H0;
    constitutiveOutputMap[NuTo::Constitutive::Output::D_INTERNAL_GRADIENT_WV_D_WV_NN_H0]                  = &evaluateData.mInternalGradientWV_dWV_NN_H0;

    constitutiveOutputMap[NuTo::Constitutive::Output::D_INTERNAL_GRADIENT_RH_D_RH_NN_H1]                  = &evaluateData.mInternalGradientRH_dRH_NN_H1;
    constitutiveOutputMap[NuTo::Constitutive::Output::D_INTERNAL_GRADIENT_RH_D_WV_NN_H1]                  = &evaluateData.mInternalGradientRH_dWV_NN_H1;
    constitutiveOutputMap[NuTo::Constitutive::Output::D_INTERNAL_GRADIENT_WV_D_WV_NN_H1]                  = &evaluateData.mInternalGradientWV_dWV_NN_H1;
    return constitutiveOutputMap;
}


bool CheckResultInRelTolerance(std::string rResultName, double rCalcValue, double rTgtValue, double rTol = 1e-6)
{
    //VHIRTHAMTODO Output formattierung für double setzen und wieder rückgängig
    double quotient = rCalcValue / rTgtValue;
    if( quotient < 1.0-rTol || quotient > 1.0+rTol)
    {
        std::cout << rResultName << " - Target Value: " << rTgtValue << " - Calculated Value: "<< rCalcValue << " - WRONG"<< std::endl;
        return false;
    }
    else
    {
        std::cout << rResultName << " - Target Value: " << rTgtValue << " - Calculated Value: "<< rCalcValue << " - CORRECT"<< std::endl;
        return true;
    }
}

template <int TDim>
void ConstitutiveOutputTest_CheckAndPrintResults(const ConstitutiveOutputTest_EvaluateData<TDim>& evaluateData)
{
 bool ValuesCorrect = true;

 // Internal Gradient
 // %%%%%%%%%%%%%%%%%
 for(unsigned int i=0; i<TDim; ++i)
 {
     ValuesCorrect = CheckResultInRelTolerance("InternalGradientRH_B ["+ std::to_string(i) +"]",evaluateData.mInternalGradientRH_B[i],-6e-13) && ValuesCorrect;
 }
 ValuesCorrect = CheckResultInRelTolerance("InternalGradientRH_N",evaluateData.mInternalGradientRH_N[0],-0.03500005) && ValuesCorrect;
 for(unsigned int i=0; i<TDim; ++i)
 {
     ValuesCorrect = CheckResultInRelTolerance("InternalGradientWV_B ["+ std::to_string(i) +"]",evaluateData.mInternalGradientWV_B[i],3.2e-7) && ValuesCorrect;
 }
 ValuesCorrect = CheckResultInRelTolerance("InternalGradientWV_N",evaluateData.mInternalGradientWV_N[0],100.00000005) && ValuesCorrect;

 // Hessian 0
 // %%%%%%%%%
 ValuesCorrect = CheckResultInRelTolerance("mInternalGradientRH_dRH_BB_H0",evaluateData.mInternalGradientRH_dRH_BB_H0[0],-3e-12) && ValuesCorrect;
 ValuesCorrect = CheckResultInRelTolerance("mInternalGradientRH_dRH_NN_H0",evaluateData.mInternalGradientRH_dRH_NN_H0[0],-0.01999995) && ValuesCorrect;
 for(unsigned int i=0; i<TDim; ++i)
 {
     ValuesCorrect = CheckResultInRelTolerance("mInternalGradientRH_dWV_BN_H0 ["+ std::to_string(i) +"]",evaluateData.mInternalGradientRH_dWV_BN_H0[i],8e-13) && ValuesCorrect;
 }
 ValuesCorrect = CheckResultInRelTolerance("mInternalGradientRH_dWV_NN_H0",evaluateData.mInternalGradientRH_dWV_NN_H0[0],0.0199999) && ValuesCorrect;
 ValuesCorrect = CheckResultInRelTolerance("mInternalGradientWV_dRH_NN_H0",evaluateData.mInternalGradientWV_dRH_NN_H0[0],-5e-8) && ValuesCorrect;
 ValuesCorrect = CheckResultInRelTolerance("mInternalGradientWV_dWV_BB_H0",evaluateData.mInternalGradientWV_dWV_BB_H0[0],1.6e-6) && ValuesCorrect;
 for(unsigned int i=0; i<TDim; ++i)
 {
     ValuesCorrect = CheckResultInRelTolerance("mInternalGradientWV_dWV_BN_H0 ["+ std::to_string(i) +"]",evaluateData.mInternalGradientWV_dWV_BN_H0[i],6.4e-7) && ValuesCorrect;
 }
 ValuesCorrect = CheckResultInRelTolerance("mInternalGradientWV_dWV_NN_H0",evaluateData.mInternalGradientWV_dWV_NN_H0[0],1e-7) && ValuesCorrect;

 // Hessian 1
 // %%%%%%%%%
  ValuesCorrect = CheckResultInRelTolerance("mInternalGradientRH_dRH_NN_H1",evaluateData.mInternalGradientRH_dRH_NN_H1[0],-0.15) && ValuesCorrect;
  ValuesCorrect = CheckResultInRelTolerance("mInternalGradientRH_dWV_NN_H1",evaluateData.mInternalGradientRH_dWV_NN_H1[0],-0.2) && ValuesCorrect;
  ValuesCorrect = CheckResultInRelTolerance("mInternalGradientWV_dWV_NN_H1",evaluateData.mInternalGradientWV_dWV_NN_H1[0],1000.0) && ValuesCorrect;

 if (!ValuesCorrect)
 {
    throw NuTo::Exception(__PRETTY_FUNCTION__,"One or more constitutive output values not correct");
 }
}

template<int TDim>
NuTo::ConstitutiveInputMap ConstitutiveOutputTest_CreateConstitutiveInputMapBoundary(ConstitutiveOutputTest_EvaluateData<TDim>& evaluateData)
{
    NuTo::ConstitutiveInputMap constitutiveInputMap;
    constitutiveInputMap[NuTo::Constitutive::Input::RELATIVE_HUMIDITY]              = &(evaluateData.mRelativeHumidity);
    constitutiveInputMap[NuTo::Constitutive::Input::WATER_VOLUME_FRACTION]          = &(evaluateData.mWaterVolumeFraction);
    return constitutiveInputMap;
}


template<int TDim>
NuTo::ConstitutiveOutputMap ConstitutiveOutputTest_CreateConstitutiveOutputMapBoundary(ConstitutiveOutputTest_EvaluateData<TDim>& evaluateData)
{
    NuTo::ConstitutiveOutputMap constitutiveOutputMap;
    constitutiveOutputMap[NuTo::Constitutive::Output::INTERNAL_GRADIENT_RELATIVE_HUMIDITY_BOUNDARY_N]     = &evaluateData.mInternalGradientRH_Boundary_N;
    constitutiveOutputMap[NuTo::Constitutive::Output::INTERNAL_GRADIENT_WATER_VOLUME_FRACTION_BOUNDARY_N] = &evaluateData.mInternalGradientWV_Boundary_N;

    constitutiveOutputMap[NuTo::Constitutive::Output::D_INTERNAL_GRADIENT_RH_D_RH_BOUNDARY_NN_H0]         = &evaluateData.mInternalGradientRH_dRH_Boundary_NN_H0;
    constitutiveOutputMap[NuTo::Constitutive::Output::D_INTERNAL_GRADIENT_WV_D_WV_BOUNDARY_NN_H0]         = &evaluateData.mInternalGradientWV_dWV_Boundary_NN_H0;
    return constitutiveOutputMap;
}


template <int TDim>
void ConstitutiveOutputTest_CheckAndPrintResultsBoundary(const ConstitutiveOutputTest_EvaluateData<TDim>& evaluateData)
{
     bool ValuesCorrect = true;

     ValuesCorrect = CheckResultInRelTolerance("mInternalGradientRH_Boundary_N",evaluateData.mInternalGradientRH_Boundary_N[0],6e-4) && ValuesCorrect;
     ValuesCorrect = CheckResultInRelTolerance("mInternalGradientWV_Boundary_N",evaluateData.mInternalGradientWV_Boundary_N[0],8e-10) && ValuesCorrect;

     ValuesCorrect = CheckResultInRelTolerance("mInternalGradientRH_dRH_Boundary_NN_H0",evaluateData.mInternalGradientRH_dRH_Boundary_NN_H0[0],1e-3) && ValuesCorrect;
     ValuesCorrect = CheckResultInRelTolerance("mInternalGradientWV_dWV_Boundary_NN_H0",evaluateData.mInternalGradientWV_dWV_Boundary_NN_H0[0],1e-9) && ValuesCorrect;


     if (!ValuesCorrect)
     {
        throw NuTo::Exception(__PRETTY_FUNCTION__,"One or more constitutive output values for the BOUNDARY element not correct");
     }
}


template<int TDim>
void ConstitutiveOutputTests(std::map<NuTo::Node::eDof,NuTo::Interpolation::eTypeOrder> rDofIPTMap)
{
    NuTo::Structure S(TDim);
    MoistureTransportControl MT(S);

    std::string testName = std::string("ConstitutiveOutputTests")+std::to_string(TDim)+"D";

    std::cout << std::endl << "-------------------------------------"
              << std::endl << "Start test: "<< testName
              << std::endl << "-------------------------------------" << std::endl;

    ConstitutiveOutputTest_SetupMoistureTransport(MT);

    SetupStructure<TDim>(S,testName);
    int SEC = SetupSection<TDim>(S);
    int IPT = SetupInterpolationType<TDim>(S,rDofIPTMap);


    switch (TDim)
    {
    case 1:
        NuTo::MeshGenerator::MeshLineSegment(S,
                                             SEC,
                                             MT.ConstitutiveLawID,
                                             IPT,
                                             {1},
                                             {1});
        break;
    case 2:
        NuTo::MeshGenerator::MeshRectangularPlane(S,
                                                  SEC,
                                                  MT.ConstitutiveLawID,
                                                  IPT,
                                                  {1,1},
                                                  {1,1});
        break;
    case 3:
        NuTo::MeshGenerator::MeshCuboid(S,
                                        SEC,
                                        MT.ConstitutiveLawID,
                                        IPT,
                                        {1,1,1},
                                        {1,1,1});
        break;

    default:
            throw NuTo::Exception(__PRETTY_FUNCTION__,"invalid dimension");
    }
    S.ElementTotalConvertToInterpolationType();
    ConstitutiveOutputTest_SetNodalValues(S,MT);




    MT.SetupStaticData();
    S.NodeBuildGlobalDofs();

    NuTo::ElementBase& element = *S.ElementGetElementPtr(0);
    NuTo::ConstitutiveBase& constLaw = *element.GetConstitutiveLaw(0);

    // %%%%%%%%%%%%%%%%%%%%%%%%
    // Constitutive law outputs
    // %%%%%%%%%%%%%%%%%%%%%%%%


    ConstitutiveOutputTest_EvaluateData<TDim> evaluateData;
    NuTo::ConstitutiveInputMap  constitutiveInputMap   = ConstitutiveOutputTest_CreateConstitutiveInputMap<TDim>(evaluateData);
    NuTo::ConstitutiveOutputMap  constitutiveOutputMap  = ConstitutiveOutputTest_CreateConstitutiveOutputMap<TDim>(evaluateData);

    evaluateData.mRelativeHumidity[0]               = 1.0;
    evaluateData.mRelativeHumidity_dt1[0]           = 0.1;
    for(unsigned int i=0; i<TDim; ++i)
    {
        evaluateData.mRelativeHumidity_Gradient[i]      = 0.2;
    }
    evaluateData.mWaterVolumeFraction[0]            = 1.0;
    evaluateData.mWaterVolumeFraction_dt1[0]        = 0.1;
    for(unsigned int i=0; i<TDim; ++i)
    {
        evaluateData.mWaterVolumeFraction_Gradient[i]   = 0.2;
    }

    constLaw.Evaluate<TDim>(&element,0,constitutiveInputMap,constitutiveOutputMap);


    ConstitutiveOutputTest_CheckAndPrintResults(evaluateData);


    // %%%%%%%%%%%%
    // Test element
    // %%%%%%%%%%%%

//    std::map<NuTo::Element::eOutput, std::shared_ptr<NuTo::ElementOutputBase>> elementOutputMap;

//    std::shared_ptr<NuTo::ElementOutputBase> internalGradient = std::make_shared<NuTo::ElementOutputBlockVectorDouble>(S.GetDofStatus());
//    std::shared_ptr<NuTo::ElementOutputBase> hessian0 = std::make_shared<NuTo::ElementOutputBlockMatrixDouble>(S.GetDofStatus());
//    std::shared_ptr<NuTo::ElementOutputBase> hessian1 = std::make_shared<NuTo::ElementOutputBlockMatrixDouble>(S.GetDofStatus());

//    elementOutputMap[NuTo::Element::eOutput::INTERNAL_GRADIENT] = internalGradient;
//    elementOutputMap[NuTo::Element::eOutput::HESSIAN_0_TIME_DERIVATIVE] = hessian0;
//    elementOutputMap[NuTo::Element::eOutput::HESSIAN_1_TIME_DERIVATIVE] = hessian1;
//    element.Evaluate(NuTo::ConstitutiveInputMap(), elementOutputMap);

//    std::cout << std::endl << "internal Gradient:" << std::endl << internalGradient->GetBlockFullVectorDouble() << std::endl;
//    std::cout << std::endl << "Hessian 0:" << std::endl << hessian0->GetBlockFullMatrixDouble() << std::endl;
//    std::cout << std::endl << "Hessian 1:" << std::endl << hessian1->GetBlockFullMatrixDouble() << std::endl;

    //S.ElementCheckHessian0(1e-10,1e-10);
//    S.ElementGetElementPtr(0)

    // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    // Constitutive law outputs - boundary element
    // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    auto GetBoundaryNodesLambda = [](NuTo::NodeBase* rNodePtr) -> bool
                                {
                                    double Tol = 1.e-6;

                                    if (rNodePtr->GetNumCoordinates()>0)
                                    {
                                        double x = rNodePtr->GetCoordinate(0);

                                        if (x >= 0.0  -Tol   && x <= 0.0  + Tol)
                                        {
                                            return true;
                                        }
                                    }
                                    return false;
                                };  // GetBoundaryNodesLambda

    int nGrpBE = S.GroupCreate("NODES");
    S.GroupAddNodeFunction(nGrpBE,GetBoundaryNodesLambda);

    int eGrpBE = S.GroupCreate("ELEMENTS");
    S.GroupAddElementsFromNodes(eGrpBE, nGrpBE, false);

    int boundaryNodeID = S.NodeCreateDOFs("RELATIVEHUMIDITY");
    NuTo::NodeBase* BEPtr = S.NodeGetNodePtr(boundaryNodeID);
    BEPtr->SetRelativeHumidity(MT.BoundaryEnvironmentalRH);
//    rS.BoundaryElementsCreate(eGrpBE,nGrpBE,rS.NodeGetNodePtr(boundaryNodeID));
    S.BoundaryElementsCreate(eGrpBE,nGrpBE,BEPtr);

      NuTo::ElementBase& boundaryElement = *S.ElementGetElementPtr(1);
      NuTo::ConstitutiveBase& boundaryConstLaw = *element.GetConstitutiveLaw(0);


      NuTo::ConstitutiveInputMap  constitutiveInputMapBoundary   = ConstitutiveOutputTest_CreateConstitutiveInputMapBoundary<TDim>(evaluateData);
      NuTo::ConstitutiveOutputMap  constitutiveOutputMapBoundary = ConstitutiveOutputTest_CreateConstitutiveOutputMapBoundary<TDim>(evaluateData);

      MT.SetupStaticData(); //VHIRTHAMTODO Check static data Coeffs in contitutive law or somewhere else --- maybe assert!!!

      boundaryConstLaw.Evaluate<TDim>(&boundaryElement,0,constitutiveInputMapBoundary,constitutiveOutputMapBoundary);


      ConstitutiveOutputTest_CheckAndPrintResultsBoundary(evaluateData);

//      boundaryElement.Evaluate(NuTo::ConstitutiveInputMap(), elementOutputMap);

//      std::cout << std::endl << "internal Gradient (boundary element):" << std::endl << internalGradient->GetBlockFullVectorDouble() << std::endl;
//      std::cout << std::endl << "Hessian 0 (boundary element):" << std::endl << hessian0->GetBlockFullMatrixDouble() << std::endl;



}

