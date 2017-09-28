#pragma once


#include "base/Timer.h"
#include "base/Exception.h"
#include "mechanics/constitutive/ConstitutiveBase.h"
#include "mechanics/constitutive/inputoutput/ConstitutiveIOBase.h"
#include "mechanics/constitutive/inputoutput/ConstitutiveIOMap.h"
#include "mechanics/constitutive/inputoutput/ConstitutiveScalar.h"
#include "mechanics/constitutive/inputoutput/ConstitutiveVector.h"
#include "mechanics/constitutive/staticData/DataMoistureTransport.h"
#include "mechanics/elements/ElementOutputBlockVectorDouble.h"
#include "mechanics/elements/ElementOutputBlockMatrixDouble.h"
#include "mechanics/structures/unstructured/Structure.h"
#include "MoistureTransport_Setup.h"

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
    for (int i = 0; i < rS.GetNumNodes(); i++)
    {

        NuTo::NodeBase& node = *rS.NodeGetNodePtr(i);

        double factor = (node.Get(NuTo::Node::eDof::COORDINATES)[0] < 0) ? 0.5 : 1.0;

        if (node.GetNum(NuTo::Node::eDof::RELATIVEHUMIDITY) != 0)
        {
            node.Set(NuTo::Node::eDof::RELATIVEHUMIDITY, 0, rMT.InitialRelativeHumidity * factor);
            node.Set(NuTo::Node::eDof::RELATIVEHUMIDITY, 1, rMT.InitialRelativeHumidity * factor * 0.01);
        }
        if (node.GetNum(NuTo::Node::eDof::WATERVOLUMEFRACTION) != 0)
        {
            node.Set(NuTo::Node::eDof::WATERVOLUMEFRACTION, 0, rMT.InitialWaterVolumeFraction * factor);
            node.Set(NuTo::Node::eDof::WATERVOLUMEFRACTION, 1, rMT.InitialWaterVolumeFraction * factor * 0.01);
        }
    }
}

template <int TDim>
NuTo::ConstitutiveInputMap ConstitutiveOutputTest_CreateConstitutiveInputMap()
{
    using namespace NuTo::Constitutive;
    NuTo::ConstitutiveInputMap constitutiveInputMap;
    std::vector<NuTo::Constitutive::eInput> inputs{
            eInput::RELATIVE_HUMIDITY,     eInput::RELATIVE_HUMIDITY_DT1,     eInput::RELATIVE_HUMIDITY_GRADIENT,
            eInput::WATER_VOLUME_FRACTION, eInput::WATER_VOLUME_FRACTION_DT1, eInput::WATER_VOLUME_FRACTION_GRADIENT};

    for (auto input : inputs)
    {
        constitutiveInputMap[input] = NuTo::ConstitutiveIOBase::makeConstitutiveIO<TDim>(input);
    }

    (*constitutiveInputMap.at(eInput::RELATIVE_HUMIDITY))[0] = 1.0;
    (*constitutiveInputMap.at(eInput::RELATIVE_HUMIDITY_DT1))[0] = 0.1;
    for (unsigned int i = 0; i < TDim; ++i)
    {
        (*constitutiveInputMap.at(eInput::RELATIVE_HUMIDITY_GRADIENT))[i] = 0.2;
    }
    (*constitutiveInputMap.at(eInput::WATER_VOLUME_FRACTION))[0] = 1.0;
    (*constitutiveInputMap.at(eInput::WATER_VOLUME_FRACTION_DT1))[0] = 0.1;
    for (unsigned int i = 0; i < TDim; ++i)
    {
        (*constitutiveInputMap.at(eInput::WATER_VOLUME_FRACTION_GRADIENT))[i] = 0.2;
    }
    return constitutiveInputMap;
}

template <int TDim>
NuTo::ConstitutiveOutputMap ConstitutiveOutputTest_CreateConstitutiveOutputMap()
{
    using namespace NuTo::Constitutive;
    NuTo::ConstitutiveOutputMap constitutiveOutputMap;
    std::vector<eOutput> outputs{
            eOutput::INTERNAL_GRADIENT_RELATIVE_HUMIDITY_B,     eOutput::INTERNAL_GRADIENT_RELATIVE_HUMIDITY_N,
            eOutput::INTERNAL_GRADIENT_WATER_VOLUME_FRACTION_B, eOutput::INTERNAL_GRADIENT_WATER_VOLUME_FRACTION_N,
            eOutput::D_INTERNAL_GRADIENT_RH_D_RH_BB_H0,         eOutput::D_INTERNAL_GRADIENT_RH_D_RH_NN_H0,
            eOutput::D_INTERNAL_GRADIENT_RH_D_WV_BN_H0,         eOutput::D_INTERNAL_GRADIENT_RH_D_WV_NN_H0,
            eOutput::D_INTERNAL_GRADIENT_WV_D_RH_NN_H0,         eOutput::D_INTERNAL_GRADIENT_WV_D_WV_BB_H0,
            eOutput::D_INTERNAL_GRADIENT_WV_D_WV_BN_H0,         eOutput::D_INTERNAL_GRADIENT_WV_D_WV_NN_H0,
            eOutput::D_INTERNAL_GRADIENT_RH_D_RH_NN_H1,         eOutput::D_INTERNAL_GRADIENT_RH_D_WV_NN_H1,
            eOutput::D_INTERNAL_GRADIENT_WV_D_WV_NN_H1};

    for (auto output : outputs)
    {
        constitutiveOutputMap[output] = NuTo::ConstitutiveIOBase::makeConstitutiveIO<TDim>(output);
    }

    return constitutiveOutputMap;
}


bool CheckResultInRelTolerance(std::string rResultName, double rCalcValue, double rTgtValue, double rTol = 1e-6)
{
    // VHIRTHAMTODO Output formattierung für double setzen und wieder rückgängig
    double quotient = rCalcValue / rTgtValue;
    if (quotient < 1.0 - rTol || quotient > 1.0 + rTol)
    {
        std::cout << rResultName << " - Target Value: " << rTgtValue << " - Calculated Value: " << rCalcValue
                  << " - WRONG" << std::endl;
        return false;
    }
    else
    {
        std::cout << rResultName << " - Target Value: " << rTgtValue << " - Calculated Value: " << rCalcValue
                  << " - CORRECT" << std::endl;
        return true;
    }
}

template <int TDim>
void ConstitutiveOutputTest_CheckAndPrintResults(const NuTo::ConstitutiveOutputMap& outputMap)
{
    using namespace NuTo::Constitutive;
    bool ValuesCorrect = true;

    // Internal Gradient
    // %%%%%%%%%%%%%%%%%
    for (unsigned int i = 0; i < TDim; ++i)
    {
        ValuesCorrect =
                CheckResultInRelTolerance("InternalGradientRH_B [" + std::to_string(i) + "]",
                                          (*outputMap.at(eOutput::INTERNAL_GRADIENT_RELATIVE_HUMIDITY_B))[i], -6e-13) &&
                ValuesCorrect;
    }
    ValuesCorrect = CheckResultInRelTolerance("InternalGradientRH_N",
                                              (*outputMap.at(eOutput::INTERNAL_GRADIENT_RELATIVE_HUMIDITY_N))[0],
                                              -0.03500005) &&
                    ValuesCorrect;
    for (unsigned int i = 0; i < TDim; ++i)
    {
        ValuesCorrect = CheckResultInRelTolerance(
                                "InternalGradientWV_B [" + std::to_string(i) + "]",
                                (*outputMap.at(eOutput::INTERNAL_GRADIENT_WATER_VOLUME_FRACTION_B))[i], 3.2e-7) &&
                        ValuesCorrect;
    }
    ValuesCorrect = CheckResultInRelTolerance("InternalGradientWV_N",
                                              (*outputMap.at(eOutput::INTERNAL_GRADIENT_WATER_VOLUME_FRACTION_N))[0],
                                              100.00000005) &&
                    ValuesCorrect;

    // Hessian 0
    // %%%%%%%%%
    ValuesCorrect = CheckResultInRelTolerance("mInternalGradientRH_dRH_BB_H0",
                                              (*outputMap.at(eOutput::D_INTERNAL_GRADIENT_RH_D_RH_BB_H0))[0], -3e-12) &&
                    ValuesCorrect;
    ValuesCorrect =
            CheckResultInRelTolerance("mInternalGradientRH_dRH_NN_H0",
                                      (*outputMap.at(eOutput::D_INTERNAL_GRADIENT_RH_D_RH_NN_H0))[0], -0.01999995) &&
            ValuesCorrect;
    for (unsigned int i = 0; i < TDim; ++i)
    {
        ValuesCorrect =
                CheckResultInRelTolerance("mInternalGradientRH_dWV_BN_H0 [" + std::to_string(i) + "]",
                                          (*outputMap.at(eOutput::D_INTERNAL_GRADIENT_RH_D_WV_BN_H0))[i], 8e-13) &&
                ValuesCorrect;
    }
    ValuesCorrect =
            CheckResultInRelTolerance("mInternalGradientRH_dWV_NN_H0",
                                      (*outputMap.at(eOutput::D_INTERNAL_GRADIENT_RH_D_WV_NN_H0))[0], 0.0199999) &&
            ValuesCorrect;
    ValuesCorrect = CheckResultInRelTolerance("mInternalGradientWV_dRH_NN_H0",
                                              (*outputMap.at(eOutput::D_INTERNAL_GRADIENT_WV_D_RH_NN_H0))[0], -5e-8) &&
                    ValuesCorrect;
    ValuesCorrect = CheckResultInRelTolerance("mInternalGradientWV_dWV_BB_H0",
                                              (*outputMap.at(eOutput::D_INTERNAL_GRADIENT_WV_D_WV_BB_H0))[0], 1.6e-6) &&
                    ValuesCorrect;
    for (unsigned int i = 0; i < TDim; ++i)
    {
        ValuesCorrect =
                CheckResultInRelTolerance("mInternalGradientWV_dWV_BN_H0 [" + std::to_string(i) + "]",
                                          (*outputMap.at(eOutput::D_INTERNAL_GRADIENT_WV_D_WV_BN_H0))[i], 6.4e-7) &&
                ValuesCorrect;
    }
    ValuesCorrect = CheckResultInRelTolerance("mInternalGradientWV_dWV_NN_H0",
                                              (*outputMap.at(eOutput::D_INTERNAL_GRADIENT_WV_D_WV_NN_H0))[0], 1e-7) &&
                    ValuesCorrect;

    // Hessian 1
    // %%%%%%%%%
    ValuesCorrect = CheckResultInRelTolerance("mInternalGradientRH_dRH_NN_H1",
                                              (*outputMap.at(eOutput::D_INTERNAL_GRADIENT_RH_D_RH_NN_H1))[0], -0.15) &&
                    ValuesCorrect;
    ValuesCorrect = CheckResultInRelTolerance("mInternalGradientRH_dWV_NN_H1",
                                              (*outputMap.at(eOutput::D_INTERNAL_GRADIENT_RH_D_WV_NN_H1))[0], -0.2) &&
                    ValuesCorrect;
    ValuesCorrect = CheckResultInRelTolerance("mInternalGradientWV_dWV_NN_H1",
                                              (*outputMap.at(eOutput::D_INTERNAL_GRADIENT_WV_D_WV_NN_H1))[0], 1000.0) &&
                    ValuesCorrect;

    if (!ValuesCorrect)
    {
        throw NuTo::Exception(__PRETTY_FUNCTION__, "One or more constitutive output values not correct");
    }
}

template <int TDim>
NuTo::ConstitutiveInputMap ConstitutiveOutputTest_CreateConstitutiveInputMapBoundary()
{
    using namespace NuTo::Constitutive;
    NuTo::ConstitutiveInputMap constitutiveInputMap;
    constitutiveInputMap[eInput::RELATIVE_HUMIDITY] =
            NuTo::ConstitutiveIOBase::makeConstitutiveIO<TDim>(eInput::RELATIVE_HUMIDITY);
    constitutiveInputMap[eInput::RELATIVE_HUMIDITY_BOUNDARY] =
            NuTo::ConstitutiveIOBase::makeConstitutiveIO<TDim>(eInput::RELATIVE_HUMIDITY_BOUNDARY);
    constitutiveInputMap[eInput::WATER_VOLUME_FRACTION] =
            NuTo::ConstitutiveIOBase::makeConstitutiveIO<TDim>(eInput::WATER_VOLUME_FRACTION);

    (*constitutiveInputMap.at(eInput::RELATIVE_HUMIDITY))[0] = 1.0;
    (*constitutiveInputMap.at(eInput::RELATIVE_HUMIDITY_BOUNDARY))[0] = 0.4;
    (*constitutiveInputMap.at(eInput::WATER_VOLUME_FRACTION))[0] = 1.0;

    return constitutiveInputMap;
}


template <int TDim>
NuTo::ConstitutiveOutputMap ConstitutiveOutputTest_CreateConstitutiveOutputMapBoundary()
{
    using namespace NuTo::Constitutive;
    NuTo::ConstitutiveOutputMap constitutiveOutputMap;
    constitutiveOutputMap[eOutput::INTERNAL_GRADIENT_RELATIVE_HUMIDITY_BOUNDARY_N] =
            NuTo::ConstitutiveIOBase::makeConstitutiveIO<TDim>(eOutput::INTERNAL_GRADIENT_RELATIVE_HUMIDITY_BOUNDARY_N);
    constitutiveOutputMap[eOutput::INTERNAL_GRADIENT_WATER_VOLUME_FRACTION_BOUNDARY_N] =
            NuTo::ConstitutiveIOBase::makeConstitutiveIO<TDim>(
                    eOutput::INTERNAL_GRADIENT_WATER_VOLUME_FRACTION_BOUNDARY_N);

    constitutiveOutputMap[eOutput::D_INTERNAL_GRADIENT_RH_D_RH_BOUNDARY_NN_H0] =
            NuTo::ConstitutiveIOBase::makeConstitutiveIO<TDim>(eOutput::D_INTERNAL_GRADIENT_RH_D_RH_BOUNDARY_NN_H0);
    constitutiveOutputMap[eOutput::D_INTERNAL_GRADIENT_WV_D_WV_BOUNDARY_NN_H0] =
            NuTo::ConstitutiveIOBase::makeConstitutiveIO<TDim>(eOutput::D_INTERNAL_GRADIENT_WV_D_WV_BOUNDARY_NN_H0);
    return constitutiveOutputMap;
}


template <int TDim>
void ConstitutiveOutputTest_CheckAndPrintResultsBoundary(const NuTo::ConstitutiveOutputMap& outputMap)
{
    using namespace NuTo::Constitutive;
    bool ValuesCorrect = true;

    ValuesCorrect = CheckResultInRelTolerance(
                            "mInternalGradientRH_Boundary_N",
                            (*outputMap.at(eOutput::INTERNAL_GRADIENT_RELATIVE_HUMIDITY_BOUNDARY_N))[0], 6e-4) &&
                    ValuesCorrect;
    ValuesCorrect = CheckResultInRelTolerance(
                            "mInternalGradientWV_Boundary_N",
                            (*outputMap.at(eOutput::INTERNAL_GRADIENT_WATER_VOLUME_FRACTION_BOUNDARY_N))[0], 8e-10) &&
                    ValuesCorrect;

    ValuesCorrect =
            CheckResultInRelTolerance("mInternalGradientRH_dRH_Boundary_NN_H0",
                                      (*outputMap.at(eOutput::D_INTERNAL_GRADIENT_RH_D_RH_BOUNDARY_NN_H0))[0], 1e-3) &&
            ValuesCorrect;
    ValuesCorrect =
            CheckResultInRelTolerance("mInternalGradientWV_dWV_Boundary_NN_H0",
                                      (*outputMap.at(eOutput::D_INTERNAL_GRADIENT_WV_D_WV_BOUNDARY_NN_H0))[0], 1e-9) &&
            ValuesCorrect;


    if (!ValuesCorrect)
    {
        throw NuTo::Exception(__PRETTY_FUNCTION__,
                              "One or more constitutive output values for the BOUNDARY element not correct");
    }
}


template <int TDim>
void ConstitutiveOutputTests(std::map<NuTo::Node::eDof, NuTo::Interpolation::eTypeOrder> rDofIPTMap)
{
    NuTo::Structure S(TDim);
    MoistureTransportControl MT(S);

    std::string testName = std::string("ConstitutiveOutputTests") + std::to_string(TDim) + "D";

    std::cout << std::endl
              << "-------------------------------------" << std::endl
              << "Start test: " << testName << std::endl
              << "-------------------------------------" << std::endl;

    ConstitutiveOutputTest_SetupMoistureTransport(MT);

    SetupStructure<TDim>(S, testName);
    auto SEC = SetupSection<TDim>(S);

    std::pair<int, int> meshInfo;

    switch (TDim)
    {
    case 1:
        meshInfo = NuTo::MeshGenerator::Grid(S, {1}, {1});
        break;
    case 2:
        meshInfo = NuTo::MeshGenerator::Grid(S, {1, 1}, {1, 1});
        break;
    case 3:
        meshInfo = NuTo::MeshGenerator::Grid(S, {1, 1, 1}, {1, 1, 1});
        break;
    default:
        throw NuTo::Exception(__PRETTY_FUNCTION__, "invalid dimension");
    }


    for (auto& it : rDofIPTMap)
        S.InterpolationTypeAdd(meshInfo.second, it.first, it.second);

    S.ElementGroupSetSection(meshInfo.first, SEC);
    S.ElementGroupSetConstitutiveLaw(meshInfo.first, MT.ConstitutiveLawID);

    SetupIntegrationType<TDim>(S, meshInfo.first);

    S.ElementTotalConvertToInterpolationType();
    ConstitutiveOutputTest_SetNodalValues(S, MT);


    MT.SetupStaticData();
    S.NodeBuildGlobalDofs();

    NuTo::ElementBase& element = *S.ElementGetElementPtr(0);

    // %%%%%%%%%%%%%%%%%%%%%%%%
    // Constitutive law outputs
    // %%%%%%%%%%%%%%%%%%%%%%%%


    NuTo::ConstitutiveInputMap constitutiveInputMap = ConstitutiveOutputTest_CreateConstitutiveInputMap<TDim>();
    NuTo::ConstitutiveOutputMap constitutiveOutputMap = ConstitutiveOutputTest_CreateConstitutiveOutputMap<TDim>();

    NuTo::Constitutive::IPConstitutiveLawBase& ipLaw = element.GetIPData().GetIPConstitutiveLaw(0);
    ipLaw.Evaluate<TDim>(constitutiveInputMap, constitutiveOutputMap);

    ConstitutiveOutputTest_CheckAndPrintResults<TDim>(constitutiveOutputMap);


    // %%%%%%%%%%%%
    // Test element
    // %%%%%%%%%%%%

    //    std::map<NuTo::ElementEnum::eOutput, std::shared_ptr<NuTo::ElementOutputBase>> elementOutputMap;

    //    std::shared_ptr<NuTo::ElementOutputBase> internalGradient =
    //    std::make_shared<NuTo::ElementOutputBlockVectorDouble>(S.GetDofStatus());
    //    std::shared_ptr<NuTo::ElementOutputBase> hessian0 =
    //    std::make_shared<NuTo::ElementOutputBlockMatrixDouble>(S.GetDofStatus());
    //    std::shared_ptr<NuTo::ElementOutputBase> hessian1 =
    //    std::make_shared<NuTo::ElementOutputBlockMatrixDouble>(S.GetDofStatus());

    //    elementOutputMap[NuTo::ElementEnum::eOutput::INTERNAL_GRADIENT] = internalGradient;
    //    elementOutputMap[NuTo::ElementEnum::eOutput::HESSIAN_0_TIME_DERIVATIVE] = hessian0;
    //    elementOutputMap[NuTo::ElementEnum::eOutput::HESSIAN_1_TIME_DERIVATIVE] = hessian1;
    //    element.Evaluate(NuTo::ConstitutiveInputMap(), elementOutputMap);

    //    std::cout << std::endl << "internal Gradient:" << std::endl << internalGradient->GetBlockFullVectorDouble() <<
    //    std::endl;
    //    std::cout << std::endl << "Hessian 0:" << std::endl << hessian0->GetBlockFullMatrixDouble() << std::endl;
    //    std::cout << std::endl << "Hessian 1:" << std::endl << hessian1->GetBlockFullMatrixDouble() << std::endl;

    // S.ElementCheckHessian0(1e-10,1e-10);
    //    S.ElementGetElementPtr(0)

    // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    // Constitutive law outputs - boundary element
    // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    auto GetBoundaryNodesLambda = [](NuTo::NodeBase* rNodePtr) -> bool {
        double Tol = 1.e-6;

        if (rNodePtr->GetNum(NuTo::Node::eDof::COORDINATES) > 0)
        {
            double x = rNodePtr->Get(NuTo::Node::eDof::COORDINATES)[0];

            if (x >= 0.0 - Tol && x <= 0.0 + Tol)
            {
                return true;
            }
        }
        return false;
    }; // GetBoundaryNodesLambda

    int nGrpBE = S.GroupCreate("NODES");
    S.GroupAddNodeFunction(nGrpBE, GetBoundaryNodesLambda);

    int eGrpBE = S.GroupCreate("ELEMENTS");
    S.GroupAddElementsFromNodes(eGrpBE, nGrpBE, false);

    int boundaryNodeID = S.NodeCreateDOFs("RELATIVEHUMIDITY");
    NuTo::NodeBase* BEPtr = S.NodeGetNodePtr(boundaryNodeID);
    BEPtr->Set(NuTo::Node::eDof::RELATIVEHUMIDITY, MT.BoundaryEnvironmentalRH);
    //    rS.BoundaryElementsCreate(eGrpBE,nGrpBE,rS.NodeGetNodePtr(boundaryNodeID));
    S.BoundaryElementsCreate(eGrpBE, nGrpBE, BEPtr);

    NuTo::ElementBase& boundaryElement = *S.ElementGetElementPtr(1);

    NuTo::ConstitutiveInputMap constitutiveInputMapBoundary =
            ConstitutiveOutputTest_CreateConstitutiveInputMapBoundary<TDim>();
    NuTo::ConstitutiveOutputMap constitutiveOutputMapBoundary =
            ConstitutiveOutputTest_CreateConstitutiveOutputMapBoundary<TDim>();

    MT.SetupStaticData(); // VHIRTHAMTODO Check static data Coeffs in contitutive law or somewhere else --- maybe
                          // assert!!!

    NuTo::Constitutive::IPConstitutiveLawBase& ipLawBoundary = boundaryElement.GetIPData().GetIPConstitutiveLaw(0);
    ipLawBoundary.Evaluate<TDim>(constitutiveInputMapBoundary, constitutiveOutputMapBoundary);

    ConstitutiveOutputTest_CheckAndPrintResultsBoundary<TDim>(constitutiveOutputMapBoundary);

    //      boundaryElement.Evaluate(NuTo::ConstitutiveInputMap(), elementOutputMap);

    //      std::cout << std::endl << "internal Gradient (boundary element):" << std::endl <<
    //      internalGradient->GetBlockFullVectorDouble() << std::endl;
    //      std::cout << std::endl << "Hessian 0 (boundary element):" << std::endl <<
    //      hessian0->GetBlockFullMatrixDouble() << std::endl;
}
