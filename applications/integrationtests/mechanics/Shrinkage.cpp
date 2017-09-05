#include "base/Timer.h"
#include "mechanics/MechanicsEnums.h"
#include "mechanics/elements/ElementBase.h"
#include "mechanics/groups/Group.h"
#include "mechanics/nodes/NodeBase.h"
#include "mechanics/structures/unstructured/Structure.h"
#include "mechanics/timeIntegration/NewmarkDirect.h"
#include "mechanics/mesh/MeshGenerator.h"
#include "physics/PhysicalConstantsSI.h"
#include "physics/PhysicalEquationsSI.h"
#include <array>
#include <boost/foreach.hpp>
#include "mechanics/constitutive/laws/MoistureTransport.h"
#include "mechanics/constitutive/staticData/DataMoistureTransport.h"
#include "mechanics/constitutive/laws/AdditiveOutput.h"
#include "mechanics/constitutive/laws/AdditiveInputExplicit.h"
#include "mechanics/constitutive/staticData/IPAdditiveOutput.h"
#include "mechanics/sections/SectionPlane.h"
#include "mechanics/sections/SectionTruss.h"
#include "mechanics/constraints/ConstraintCompanion.h"
#include "mechanics/timeIntegration/postProcessing/PostProcessor.h"

#include "visualize/VisualizeEnum.h"


//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//  Setup Preprocessor
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


/*---------------------------------------------*\
|*                  DEFINES                    *|
\*---------------------------------------------*/


// --- Time integration scheme
// ---------------------------
#define RES_TOLERANCE_MECHANICS 1e-8
#define RES_TOLERANCE_MOISTURE_TRANSPORT 1e-17
#define MAX_ITERATION 40


// --- Processor/OpenMp
// --------------------

#ifdef _OPENMP
#define TESTNUM_PROC 4
#elif HAVE_PARDISO
#define TESTNUM_PROC 4
#else
#define TESTNUM_PROC 1
#endif

// --- Other
// ---------

#define TEST_TEMPERATURE 293.15
#define TEST_YOUNGSMODULUS 30.e9
#define TEST_MACROSCOPICBULKMODULUS 30.e9
#define TEST_SOLIDPHASEBULKMODULUS 7.5e9

/*---------------------------------------------*\
|*                  TYPEDEFS                   *|
\*---------------------------------------------*/

typedef boost::ptr_map<int, NuTo::NodeBase> NodeMap;


//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//  Setup structs
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

/*---------------------------------------------*\
|*                    time                     *|
\*---------------------------------------------*/

struct TimeControl
{
    double delta_t = 1.0 / 1.0 * 1.0 * 24.0 * 60.0 * 60.0;
    double t_write = 1.0 / 1.0 * 1.0 * 24.0 * 60.0 * 60.0;
    double t_final = 20.0 / 1.0 * 1.0 * 24.0 * 60.0 * 60.0;
    double BC_TransitionTime = 24.0 * 60.0 * 60.0;
};

struct MechanicsControl
{
    //! @brief constructor
    MechanicsControl(NuTo::Structure& rS, NuTo::ConstitutiveBase& rCL)
        : mS(rS)
        , mCL(rCL)
    {
        switch (mCL.GetType())
        {
        case NuTo::Constitutive::eConstitutiveType::LINEAR_ELASTIC_ENGINEERING_STRESS:
            break;
        default:
            throw NuTo::Exception(__PRETTY_FUNCTION__, "The referenced constitutive law is not a mechanics model");
        }
    }


    void SetParametersConstitutiveLaw()
    {
        // set variables

        mCL.SetParameterDouble(NuTo::Constitutive::eConstitutiveParameter::DENSITY, mDensity);
        mCL.SetParameterDouble(NuTo::Constitutive::eConstitutiveParameter::POISSONS_RATIO, mPoissonRatio);
        mCL.SetParameterDouble(NuTo::Constitutive::eConstitutiveParameter::YOUNGS_MODULUS, mYoungsModulus);
    }

    void AddConstraint(std::function<bool(NuTo::NodeBase*)> rGetNodeFunction, NuTo::eDirection rDirection,
                       std::function<double(double)> rDisplacementFunction = NuTo::Constraint::RhsConstant(0))
    {
        int GRPNodesConstraint = mS.GroupCreate("Nodes");
        mS.GroupAddNodeFunction(GRPNodesConstraint, rGetNodeFunction);

        const auto& group = *mS.GroupGetGroupPtr(GRPNodesConstraint)->AsGroupNode();

        mS.Constraints().Add(NuTo::Node::eDof::DISPLACEMENTS,
                             NuTo::Constraint::Component(group, {rDirection}, rDisplacementFunction));
    }

    NuTo::Structure& mS;
    NuTo::ConstitutiveBase& mCL;

    double mYoungsModulus = TEST_YOUNGSMODULUS;
    double mPoissonRatio = 0.0;
    double mDensity = 1.0;
};


struct MoistureTransportControl
{

    // ctor
    MoistureTransportControl(NuTo::Structure& rS, NuTo::ConstitutiveBase& rMT)
        : mS(rS)
        , mMT(rMT)
    {
        if (mMT.GetType() != NuTo::Constitutive::eConstitutiveType::MOISTURE_TRANSPORT)
            throw NuTo::Exception(__PRETTY_FUNCTION__,
                                  "The referenced constitutive law is not a moisture transport model");

        // Values fitted from figure in Johannessons paper
        AdsorptionCoeffs(0) = 0.0;
        AdsorptionCoeffs(1) = 0.19692057340725558;
        AdsorptionCoeffs(2) = -0.28253538941816925;
        AdsorptionCoeffs(3) = 0.22661481601091368;

        DesorptionCoeffs(0) = 0.0;
        DesorptionCoeffs(1) = 0.26719233184420238;
        DesorptionCoeffs(2) = -0.41030868184510738;
        DesorptionCoeffs(3) = 0.32511635000090505;
    }

public:
    void SetParametersConstitutiveLaw()
    {
        // set variables
        mMT.SetParameterBool(NuTo::Constitutive::eConstitutiveParameter::ENABLE_MODIFIED_TANGENTIAL_STIFFNESS,
                             EnableModiefiedTangentialStiffness);
        mMT.SetParameterBool(NuTo::Constitutive::eConstitutiveParameter::ENABLE_SORPTION_HYSTERESIS,
                             EnableSorptionHysteresis);

        mMT.SetParameterDouble(NuTo::Constitutive::eConstitutiveParameter::BOUNDARY_DIFFUSION_COEFFICIENT_RH,
                               BoundaryDiffusionCoefficientRH);
        mMT.SetParameterDouble(NuTo::Constitutive::eConstitutiveParameter::BOUNDARY_DIFFUSION_COEFFICIENT_WV,
                               BoundaryDiffusionCoefficientWV);
        mMT.SetParameterDouble(NuTo::Constitutive::eConstitutiveParameter::DENSITY_WATER, DensityWater);
        mMT.SetParameterDouble(NuTo::Constitutive::eConstitutiveParameter::DIFFUSION_COEFFICIENT_RH,
                               DiffusionCoefficientRH);
        mMT.SetParameterDouble(NuTo::Constitutive::eConstitutiveParameter::DIFFUSION_COEFFICIENT_WV,
                               DiffusionCoefficientWV);
        mMT.SetParameterDouble(NuTo::Constitutive::eConstitutiveParameter::DIFFUSION_EXPONENT_RH, DiffusionExponentRH);
        mMT.SetParameterDouble(NuTo::Constitutive::eConstitutiveParameter::DIFFUSION_EXPONENT_WV, DiffusionExponentWV);
        mMT.SetParameterDouble(NuTo::Constitutive::eConstitutiveParameter::GRADIENT_CORRECTION_ADSORPTION_DESORPTION,
                               GradientCorrectionAdsorpionDesorption);
        mMT.SetParameterDouble(NuTo::Constitutive::eConstitutiveParameter::GRADIENT_CORRECTION_DESORPTION_ADSORPTION,
                               GradientCorrectionDesorptionAdsorpion);
        mMT.SetParameterDouble(NuTo::Constitutive::eConstitutiveParameter::MASS_EXCHANGE_RATE, MassExchangeRate);
        mMT.SetParameterDouble(NuTo::Constitutive::eConstitutiveParameter::PORE_VOLUME_FRACTION, PoreVolumeFraction);
        mMT.SetParameterDouble(NuTo::Constitutive::eConstitutiveParameter::DENSITY_SATURATED_WATER_VAPOR,
                               DensitySaturatedWaterVapor);

        mMT.SetParameterFullVectorDouble(NuTo::Constitutive::eConstitutiveParameter::POLYNOMIAL_COEFFICIENTS_ADSORPTION,
                                         AdsorptionCoeffs);
        mMT.SetParameterFullVectorDouble(NuTo::Constitutive::eConstitutiveParameter::POLYNOMIAL_COEFFICIENTS_DESORPTION,
                                         DesorptionCoeffs);


        // Calculate equilibrium water volume fraction
        InitialWaterVolumeFraction = mMT.GetEquilibriumWaterVolumeFraction(InitialRelativeHumidity, DesorptionCoeffs);
    }


public:
    void ApplyInitialNodalValues()
    {
        unsigned int NNodes = mS.GetNumNodes();

        for (unsigned int i = 0; i < NNodes; i++)
        {


            if (mS.NodeGetNodePtr(i)->GetNum(NuTo::Node::eDof::RELATIVEHUMIDITY) != 0)
            {
                mS.NodeGetNodePtr(i)->Set(NuTo::Node::eDof::RELATIVEHUMIDITY, 0, InitialRelativeHumidity);
            }
            if (mS.NodeGetNodePtr(i)->GetNum(NuTo::Node::eDof::WATERVOLUMEFRACTION) != 0)
            {
                mS.NodeGetNodePtr(i)->Set(NuTo::Node::eDof::WATERVOLUMEFRACTION, 0, InitialWaterVolumeFraction);
            }
        }
    }

    void SetupStaticData()
    {
        for (int i = 0; i < mS.GetNumElements(); i++)
        {
            for (int theIP = 0; theIP < mS.ElementGetElementPtr(i)->GetNumIntegrationPoints(); theIP++)
            {
                NuTo::Constitutive::IPAdditiveOutput& ipLawAO = dynamic_cast<NuTo::Constitutive::IPAdditiveOutput&>(
                        mS.ElementGetElementPtr(i)->GetIPData().GetIPConstitutiveLaw(theIP));

                NuTo::Constitutive::StaticData::DataMoistureTransport& moistureData =
                        ipLawAO.GetSublawData<NuTo::MoistureTransport>(&mMT).GetData(); // finally the data.

                moistureData.SetLastSorptionCoeff(mMT.GetParameterFullVectorDouble(
                        NuTo::Constitutive::eConstitutiveParameter::POLYNOMIAL_COEFFICIENTS_DESORPTION));
                moistureData.SetCurrentSorptionCoeff(mMT.GetParameterFullVectorDouble(
                        NuTo::Constitutive::eConstitutiveParameter::POLYNOMIAL_COEFFICIENTS_DESORPTION));
                moistureData.SetLastRelHumValue(InitialRelativeHumidity);
                moistureData.SetDesorption(SorptionHistoryDesorption);
            }
        }
    }


    // general
    bool EnableModiefiedTangentialStiffness = false;

    double InitialRelativeHumidity = 1.0;
    double InitialWaterVolumeFraction = 0.; //! --> Calculated from relative humidity
    double MassExchangeRate = 3.42e-7;
    double PoreVolumeFraction = 0.25;
    double DiffusionCoefficientRH = 3.9e-12;
    double DiffusionExponentRH = 1.0;
    double DensitySaturatedWaterVapor = 0.0173;
    double DensityWater = 999.97;
    double DiffusionCoefficientWV = 1.17e-7;
    double DiffusionExponentWV = 2.0;

    // boundary condition
    double BoundaryEnvironmentalRH = 0.45;
    double BoundaryDiffusionCoefficientRH = 1.0e-10 * 1000;
    double BoundaryDiffusionCoefficientWV = 1.0e-7 * 1000;

    // sorption isotherms
    bool EnableSorptionHysteresis = false;
    bool SorptionHistoryDesorption = true;
    double GradientCorrectionDesorptionAdsorpion = 0.26;
    double GradientCorrectionAdsorpionDesorption = 0.56;
    Eigen::Vector4d AdsorptionCoeffs = Eigen::Vector4d::Zero();
    Eigen::Vector4d DesorptionCoeffs = Eigen::Vector4d::Zero();


    // references
    NuTo::Structure& mS;
    NuTo::ConstitutiveBase& mMT;
};


//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//  Setup Functions
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


/*---------------------------------------------*\
|*              boundary elements              *|
\*---------------------------------------------*/

template <int TDim>
void SetupConstrainedNodeBoundaryElements(NuTo::Structure& rS,
                                          std::function<bool(NuTo::NodeBase*)> rFunctionGetBoundaryNode,
                                          std::function<double(double)> rBoundaryConstraintFunction)
{
    int nGrpBE = rS.GroupCreate("NODES");
    rS.GroupAddNodeFunction(nGrpBE, rFunctionGetBoundaryNode);

    int eGrpBE = rS.GroupCreate("ELEMENTS");
    rS.GroupAddElementsFromNodes(eGrpBE, nGrpBE, false);

    std::set<NuTo::Node::eDof> controlNodeDofs;
    controlNodeDofs.insert(NuTo::Node::eDof::RELATIVEHUMIDITY);

    int boundaryControlNodeID = rS.NodeCreateDOFs(controlNodeDofs);
    NuTo::NodeBase* controlNodePtr = rS.NodeGetNodePtr(boundaryControlNodeID);
    int groupBoundaryElements = rS.BoundaryElementsCreate(eGrpBE, nGrpBE, controlNodePtr);

    rS.Constraints().Add(NuTo::Node::eDof::RELATIVEHUMIDITY,
                         NuTo::Constraint::Value(*controlNodePtr, rBoundaryConstraintFunction));

    // Set Integration type - default not sufficient
    std::vector<int> boundaryElementIDs;
    rS.ElementGroupGetMembers(groupBoundaryElements, boundaryElementIDs);


    for (int elementId : boundaryElementIDs)
    {
        NuTo::ElementBase* elementPtr = rS.ElementGetElementPtr(elementId);
        switch (TDim)
        {
        case 1:
            elementPtr->SetIntegrationType(
                    *rS.GetPtrIntegrationType(NuTo::eIntegrationType::IntegrationType0DBoundary));
            break;
        case 2:
            elementPtr->SetIntegrationType(
                    *rS.GetPtrIntegrationType(NuTo::eIntegrationType::IntegrationType1D2NGauss2Ip));
            break;
        case 3:
            elementPtr->SetIntegrationType(
                    *rS.GetPtrIntegrationType(NuTo::eIntegrationType::IntegrationType2D4NGauss4Ip));
            break;
        default:
            throw NuTo::Exception(__PRETTY_FUNCTION__, "Invalid dimension");
        }
    }
}


/*---------------------------------------------*\
|*              integration type               *|
\*---------------------------------------------*/

template <int TDim>
void SetupIntegrationType(NuTo::Structure& rS, int rIPT)
{
    switch (TDim)
    {
    case 1:
        rS.InterpolationTypeSetIntegrationType(rIPT, NuTo::eIntegrationType::IntegrationType1D2NGauss2Ip);
        break;
    case 2:
        rS.InterpolationTypeSetIntegrationType(rIPT, NuTo::eIntegrationType::IntegrationType2D4NGauss4Ip);
        break;
    case 3:
        rS.InterpolationTypeSetIntegrationType(rIPT, NuTo::eIntegrationType::IntegrationType3D8NGauss2x2x2Ip);
        break;
    default:
        throw NuTo::Exception(__PRETTY_FUNCTION__, "Invalid dimension");
    }
}


/*---------------------------------------------*\
|*            multi processor setup            *|
\*---------------------------------------------*/


inline void SetupMultiProcessor(NuTo::Structure& rS)
{

    rS.SetNumProcessors(TESTNUM_PROC);
#ifdef _OPENMP
    std::cout << "OpenMP enabled" << std::endl;
    rS.CalculateMaximumIndependentSets();
#endif
}

/*---------------------------------------------*\
|*                 section                     *|
\*---------------------------------------------*/

template <int TDim>
std::shared_ptr<NuTo::Section> SetupSection(NuTo::Structure& rS, double rAreaThickness = 1.0)
{
    switch (TDim)
    {
    case 1:
    {
        auto Sec = NuTo::SectionTruss::Create(rAreaThickness);
        rS.ElementTotalSetSection(Sec);
        return Sec;
    }

    case 2:
    {
        auto Sec = NuTo::SectionPlane::Create(rAreaThickness, false);
        rS.ElementTotalSetSection(Sec);
        return Sec;
    }

    case 3:
    {
        // there is no need to attach a section to 3D elements
        // to make this function work with arbitrary dimensions, we just attach a dummy truss
        auto Sec = NuTo::SectionTruss::Create(-42.0);
        rS.ElementTotalSetSection(Sec);
        return Sec;
    }

    default:
        throw NuTo::Exception(__PRETTY_FUNCTION__, "Invalid dimension");
    }
}

/*---------------------------------------------*\
|*                 structure                   *|
\*---------------------------------------------*/

void SetupStructure(NuTo::Structure& rS, std::string rTestName)
{
    rS.SetNumTimeDerivatives(1);
    rS.SetShowTime(false);

    NuTo::Logger& Log = rS.GetLogger();
    Log.SetQuiet(false);
    Log.OpenFile(rTestName + ".log");
}


/*---------------------------------------------*\
|*              time integration               *|
\*---------------------------------------------*/

inline void SetupTimeIntegration(NuTo::NewmarkDirect& rTI, const TimeControl& rTC, const std::string& rResultDir,
                                 bool rStaggered)
{
    rTI.SetPerformLineSearch(false);
    rTI.SetVerboseLevel(0);
    rTI.SetToleranceResidual(NuTo::Node::eDof::DISPLACEMENTS, RES_TOLERANCE_MECHANICS);
    rTI.SetToleranceResidual(NuTo::Node::eDof::RELATIVEHUMIDITY, RES_TOLERANCE_MOISTURE_TRANSPORT);
    rTI.SetToleranceResidual(NuTo::Node::eDof::WATERVOLUMEFRACTION, RES_TOLERANCE_MOISTURE_TRANSPORT);
    rTI.SetMaxNumIterations(MAX_ITERATION);

    rTI.SetTimeStep(rTC.delta_t);
    rTI.PostProcessing().SetMinTimeStepPlot(rTC.t_write);

    rTI.PostProcessing().SetResultDirectory(rResultDir, true);

    if (rStaggered)
    {
        rTI.AddCalculationStep({NuTo::Node::eDof::RELATIVEHUMIDITY, NuTo::Node::eDof::WATERVOLUMEFRACTION});
        rTI.AddCalculationStep({NuTo::Node::eDof::DISPLACEMENTS});
    }
}

/*---------------------------------------------*\
|*                 visualize                   *|
\*---------------------------------------------*/

inline void SetupVisualize(NuTo::Structure& rS, bool rVisualizeShrinkageStrains = false)
{
#ifdef ENABLE_VISUALIZE
    int visGrp = rS.GroupCreate(NuTo::eGroupId::Elements);
    rS.GroupAddElementsTotal(visGrp);
    rS.AddVisualizationComponent(visGrp, NuTo::eVisualizeWhat::DISPLACEMENTS);
    rS.AddVisualizationComponent(visGrp, NuTo::eVisualizeWhat::RELATIVE_HUMIDITY);
    rS.AddVisualizationComponent(visGrp, NuTo::eVisualizeWhat::WATER_VOLUME_FRACTION);
    rS.AddVisualizationComponent(visGrp, NuTo::eVisualizeWhat::ENGINEERING_STRAIN);
    rS.AddVisualizationComponent(visGrp, NuTo::eVisualizeWhat::ENGINEERING_STRESS);
    rS.AddVisualizationComponent(visGrp, NuTo::eVisualizeWhat::PRINCIPAL_ENGINEERING_STRESS);
    if (rVisualizeShrinkageStrains)
        rS.AddVisualizationComponent(visGrp, NuTo::eVisualizeWhat::SHRINKAGE_STRAIN);
#endif // ENABLE_VISUALIZE
}


//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//  Check Results
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

template <int TDim>
void CheckMechanicsResultsStressBased(NuTo::Structure& rS)
{
    const double temperature = TEST_TEMPERATURE;
    const double capStressFactor = NuTo::SI::DensityLiquidWater(temperature) * NuTo::SI::IdealGasConstant *
                                   temperature / NuTo::SI::MolarMassWater;

    const double capStress = 0.062035 * std::log(0.4) * capStressFactor;

    const NodeMap& nodePtrMap = rS.NodeGetNodeMap();
    for (auto it : nodePtrMap)
    {
        const NuTo::NodeBase* nodePtr = it.second;
        if (nodePtr->GetNum(NuTo::Node::eDof::DISPLACEMENTS) < 1)
        {
            continue; // Nodes without Displacements cant be checked
        }
        for (int i = 0; i < TDim; ++i)
        {
            double coord = nodePtr->Get(NuTo::Node::eDof::COORDINATES)[i];
            double refDisp = capStress / TEST_YOUNGSMODULUS * coord;
            double disp = nodePtr->Get(NuTo::Node::eDof::DISPLACEMENTS)[i];
            double diff = std::abs(refDisp) - std::abs(disp);
            const double tolerance = 1e-10;
            if ((diff > tolerance || diff < -tolerance)) // && coordX > 0)
            {
                throw NuTo::Exception(__PRETTY_FUNCTION__, "One ore more calculated Displacements are not correct");
            }
        }
    }
    std::cout << "Displacements correct!" << std::endl;
}

template <int TDim>
void CheckMoistureTransportResults(NuTo::Structure& rS, std::vector<int>, std::vector<double>)
{


    constexpr double tolerance = 0.00001; // Tolerance because not all necessary value (sorption curve) are given in the
    // paper and must be approximated
    unsigned int numMismatchingValues = 0;


    const NodeMap& nodePtrMap = rS.NodeGetNodeMap();
    for (auto it : nodePtrMap)
    {
        const NuTo::NodeBase* nodePtr = it.second;


        if (nodePtr->GetNum(NuTo::Node::eDof::WATERVOLUMEFRACTION) > 0)
        {
            double nodalWVF = nodePtr->Get(NuTo::Node::eDof::WATERVOLUMEFRACTION)[0];
            double eqWVF = 0.062035;
            if (nodalWVF < eqWVF - tolerance || nodalWVF > eqWVF + tolerance)
            {
                ++numMismatchingValues;
            }
        }
        if (nodePtr->GetNum(NuTo::Node::eDof::WATERVOLUMEFRACTION) > 0)
        {
            double nodalRH = nodePtr->Get(NuTo::Node::eDof::RELATIVEHUMIDITY)[0];
            double eqRH = 0.4;
            if (nodalRH < eqRH - tolerance || nodalRH > eqRH + tolerance)
            {
                ++numMismatchingValues;
            }
        }
    }
    if (numMismatchingValues > 0)
    {
        throw NuTo::Exception(__PRETTY_FUNCTION__, "One ore more calculated relative humidity/water volue fraction "
                                                   "values exceeds the tolerance when compared to reference values");
    }
    std::cout << "Water volume fraction and relative humidity correct!" << std::endl;
}


//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//  Tests
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


//! @brief performs a simulation in the desired dimension
//! @param rN: array with number of elements in each direction
//! @param rL: array with length of elements in each direction
template <int TDim>
void ShrinkageTestStressBased(std::vector<int> rN, std::vector<double> rL,
                              std::map<NuTo::Node::eDof, NuTo::Interpolation::eTypeOrder> rDofIPTMap,
                              bool rStaggered = false)
{

    std::string testName = std::string("StressBased") + std::to_string(TDim) + "D";
    if (rStaggered)
        testName += "_staggered";
    std::string resultDir = std::string("./Shrinkage_") + testName;

    std::cout << std::endl
              << "--------------------------------------------------------------------------" << std::endl
              << "Start test: " << testName << std::endl
              << "--------------------------------------------------------------------------" << std::endl;

    // Allocate neccessary stuff
    NuTo::Structure S(TDim);
    NuTo::NewmarkDirect TI(&S);
    int CL_LE_ID = S.ConstitutiveLawCreate(NuTo::Constitutive::eConstitutiveType::LINEAR_ELASTIC_ENGINEERING_STRESS);
    int CL_SCSB_ID = S.ConstitutiveLawCreate(NuTo::Constitutive::eConstitutiveType::SHRINKAGE_CAPILLARY_STRESS_BASED);
    int CL_MT_ID = S.ConstitutiveLawCreate(NuTo::Constitutive::eConstitutiveType::MOISTURE_TRANSPORT);
    int CL_AL_ID = S.ConstitutiveLawCreate(NuTo::Constitutive::eConstitutiveType::ADDITIVE_OUTPUT);


    NuTo::ConstitutiveBase* CL_LE_Ptr = S.ConstitutiveLawGetConstitutiveLawPtr(CL_LE_ID);
    NuTo::ConstitutiveBase* CL_SCSB_Ptr = S.ConstitutiveLawGetConstitutiveLawPtr(CL_SCSB_ID);
    NuTo::ConstitutiveBase* CL_MT_Ptr = S.ConstitutiveLawGetConstitutiveLawPtr(CL_MT_ID);
    NuTo::AdditiveOutput* CL_AL_Ptr =
            static_cast<NuTo::AdditiveOutput*>(S.ConstitutiveLawGetConstitutiveLawPtr(CL_AL_ID));

    TimeControl tCtrl;
    tCtrl.t_final = 365.0 * 24.0 * 60.0 * 60.0;
    tCtrl.delta_t = tCtrl.t_final;
    tCtrl.t_write = tCtrl.t_final;


    MoistureTransportControl MTCtrl(S, *CL_MT_Ptr);
    MTCtrl.InitialRelativeHumidity = 0.95;
    MTCtrl.BoundaryEnvironmentalRH = 0.40;
    MTCtrl.MassExchangeRate = 1.;
    MTCtrl.DiffusionCoefficientRH = 1e-5;
    MTCtrl.DiffusionCoefficientWV = 1.;
    MTCtrl.BoundaryDiffusionCoefficientRH = 1.0e-4;
    MTCtrl.BoundaryDiffusionCoefficientWV = 1.0;
    MTCtrl.SetParametersConstitutiveLaw();

    MechanicsControl MeCtrl(S, *CL_LE_Ptr);
    MeCtrl.SetParametersConstitutiveLaw();

    CL_SCSB_Ptr->SetParameterDouble(NuTo::Constitutive::eConstitutiveParameter::TEMPERATURE, TEST_TEMPERATURE);


    CL_AL_Ptr->AddConstitutiveLaw(*CL_MT_Ptr);
    CL_AL_Ptr->AddConstitutiveLaw(*CL_LE_Ptr);
    CL_AL_Ptr->AddConstitutiveLaw(*CL_SCSB_Ptr);


    SetupStructure(S, testName);
    auto SEC = SetupSection<TDim>(S);

    auto meshInfo = NuTo::MeshGenerator::Grid(S, rL, rN);

    for (auto& it : rDofIPTMap)
        S.InterpolationTypeAdd(meshInfo.second, it.first, it.second);

    S.ElementGroupSetSection(meshInfo.first, SEC);
    S.ElementGroupSetConstitutiveLaw(meshInfo.first, CL_AL_ID);

    SetupIntegrationType<TDim>(S, meshInfo.first);


    S.ElementTotalConvertToInterpolationType(); // old used values 1.0e-12,0.001
    MTCtrl.ApplyInitialNodalValues();


    auto LambdaGetBoundaryNodes = [rL](NuTo::NodeBase* rNodePtr) -> bool {
        double Tol = 1.e-6;
        if (rNodePtr->GetNum(NuTo::Node::eDof::COORDINATES) > 0)
        {
            double x = rNodePtr->Get(NuTo::Node::eDof::COORDINATES)[0];
            if ((x >= 0.0 - Tol && x <= 0.0 + Tol) || (x >= rL[0] - Tol && x <= rL[0] + Tol))
            {
                return true;
            }

            if (TDim > 1)
            {
                double y = rNodePtr->Get(NuTo::Node::eDof::COORDINATES)[1];
                if ((y >= 0.0 - Tol && y <= 0.0 + Tol) || (y >= rL[1] - Tol && y <= rL[1] + Tol))
                {
                    return true;
                }
            }
            if (TDim > 2)
            {
                double z = rNodePtr->Get(NuTo::Node::eDof::COORDINATES)[2];
                if ((z >= 0.0 - Tol && z <= 0.0 + Tol) || (z >= rL[2] - Tol && z <= rL[2] + Tol))
                {
                    return true;
                }
            }
        }
        return false;
    }; // GetBoundaryNodesLambda


    auto LambdaTimeDepBoundaryRH = [&MTCtrl, &tCtrl](double rTime) -> double {
        if (rTime == 0.0)
        {
            return MTCtrl.InitialRelativeHumidity;
        }
        else
        {
            if (rTime < tCtrl.BC_TransitionTime)
            {
                return MTCtrl.InitialRelativeHumidity -
                       sin(rTime / tCtrl.BC_TransitionTime * 3.14 / 2.0) *
                               (MTCtrl.InitialRelativeHumidity - MTCtrl.BoundaryEnvironmentalRH);
            }
            {
                return MTCtrl.BoundaryEnvironmentalRH;
            }
        }
    }; // TimeDepBoundaryRHLambda


    SetupConstrainedNodeBoundaryElements<TDim>(S, LambdaGetBoundaryNodes, LambdaTimeDepBoundaryRH);


    auto lambdaGetNodeLeftBottomFront = [rL](NuTo::NodeBase* rNodePtr) -> bool {
        if (rNodePtr->GetNum(NuTo::Node::eDof::DISPLACEMENTS) == 0)
            return false;
        double Tol = 1.e-6;
        if (rNodePtr->GetNum(NuTo::Node::eDof::COORDINATES) > 0)
        {
            double x = 0.0, y = 0.0, z = 0.0;
            x = rNodePtr->Get(NuTo::Node::eDof::COORDINATES)[0];
            if (TDim > 1)
                y = rNodePtr->Get(NuTo::Node::eDof::COORDINATES)[1];
            if (TDim > 2)
                z = rNodePtr->Get(NuTo::Node::eDof::COORDINATES)[2];

            if (x >= 0.0 - Tol && x <= 0.0 + Tol && y >= 0.0 - Tol && y <= 0.0 + Tol && z >= 0.0 - Tol &&
                z <= 0.0 + Tol)
            {
                return true;
            }
        }
        return false;
    }; // lambdaGetNodeLeftBottom

    auto lambdaGetNodeLeftTopFront = [rL](NuTo::NodeBase* rNodePtr) -> bool {
        if (rNodePtr->GetNum(NuTo::Node::eDof::DISPLACEMENTS) == 0)
            return false;
        double Tol = 1.e-6;
        if (rNodePtr->GetNum(NuTo::Node::eDof::COORDINATES) > 0)
        {
            double x = 0.0, y = 0.0, z = 0.0;
            x = rNodePtr->Get(NuTo::Node::eDof::COORDINATES)[0];
            if (TDim > 1)
                y = rNodePtr->Get(NuTo::Node::eDof::COORDINATES)[1];
            if (TDim > 2)
                z = rNodePtr->Get(NuTo::Node::eDof::COORDINATES)[2];

            if (x >= 0.0 - Tol && x <= 0.0 + Tol && y >= 0.0 - Tol && y <= 0.0 + Tol && z >= rL[2] - Tol &&
                z <= rL[2] + Tol)
            {
                return true;
            }
        }
        return false;
    }; // lambdaGetNodeLeftBottom

    auto lambdaGetNodeLeftBottomBack = [rL](NuTo::NodeBase* rNodePtr) -> bool {
        if (rNodePtr->GetNum(NuTo::Node::eDof::DISPLACEMENTS) == 0)
            return false;
        double Tol = 1.e-6;
        if (rNodePtr->GetNum(NuTo::Node::eDof::COORDINATES) > 0)
        {
            double x = 0.0, y = 0.0, z = 0.0;
            x = rNodePtr->Get(NuTo::Node::eDof::COORDINATES)[0];
            if (TDim > 1)
                y = rNodePtr->Get(NuTo::Node::eDof::COORDINATES)[1];
            if (TDim > 2)
                z = rNodePtr->Get(NuTo::Node::eDof::COORDINATES)[2];

            if (x >= 0.0 - Tol && x <= 0.0 + Tol && y >= rL[1] - Tol && y <= rL[1] + Tol && z >= 0.0 - Tol &&
                z <= 0.0 + Tol)
            {
                return true;
            }
        }
        return false;
    }; // lambdaGetNodeLeftBottom


    MeCtrl.AddConstraint(lambdaGetNodeLeftBottomFront, NuTo::eDirection::X);
    if (TDim > 1)
    {
        MeCtrl.AddConstraint(lambdaGetNodeLeftBottomFront, NuTo::eDirection::Y);
        MeCtrl.AddConstraint(lambdaGetNodeLeftBottomBack, NuTo::eDirection::X);
    }

    if (TDim > 2)
    {
        MeCtrl.AddConstraint(lambdaGetNodeLeftBottomFront, NuTo::eDirection::Z);
        MeCtrl.AddConstraint(lambdaGetNodeLeftBottomBack, NuTo::eDirection::Z);
        MeCtrl.AddConstraint(lambdaGetNodeLeftTopFront, NuTo::eDirection::X);
        MeCtrl.AddConstraint(lambdaGetNodeLeftTopFront, NuTo::eDirection::Y);
    }


    MTCtrl.SetupStaticData();
    S.NodeBuildGlobalDofs();

    SetupMultiProcessor(S);

    SetupVisualize(S);

    SetupTimeIntegration(TI, tCtrl, resultDir, rStaggered);
    NuTo::Timer timer("shrinkagetest");
    TI.Solve(tCtrl.t_final);

    CheckMoistureTransportResults<TDim>(S, rN, rL);
    CheckMechanicsResultsStressBased<TDim>(S);
}


//! @brief performs a simulation in the desired dimension
//! @param rN: array with number of elements in each direction
//! @param rL: array with length of elements in each direction
template <int TDim>
void ShrinkageTestStrainBased(std::vector<int> rN, std::vector<double> rL,
                              std::map<NuTo::Node::eDof, NuTo::Interpolation::eTypeOrder> rDofIPTMap,
                              bool rStaggered = false)
{
    std::string testName = std::string("StrainBased") + std::to_string(TDim) + "D";
    if (rStaggered)
        testName += "_staggered";
    std::string resultDir = std::string("./Shrinkage_") + testName;

    std::cout << std::endl
              << "--------------------------------------------------------------------------" << std::endl
              << "Start test: " << testName << std::endl
              << "--------------------------------------------------------------------------" << std::endl;

    // Allocate neccessary stuff
    NuTo::Structure S(TDim);
    NuTo::NewmarkDirect TI(&S);
    int CL_LE_ID = S.ConstitutiveLawCreate(NuTo::Constitutive::eConstitutiveType::LINEAR_ELASTIC_ENGINEERING_STRESS);
    int CL_SCSB_ID = S.ConstitutiveLawCreate(NuTo::Constitutive::eConstitutiveType::SHRINKAGE_CAPILLARY_STRAIN_BASED);
    int CL_AIE_ID = S.ConstitutiveLawCreate(NuTo::Constitutive::eConstitutiveType::ADDITIVE_INPUT_EXPLICIT);
    int CL_MT_ID = S.ConstitutiveLawCreate(NuTo::Constitutive::eConstitutiveType::MOISTURE_TRANSPORT);
    int CL_AO_ID = S.ConstitutiveLawCreate(NuTo::Constitutive::eConstitutiveType::ADDITIVE_OUTPUT);


    auto CL_LE_Ptr = S.ConstitutiveLawGetConstitutiveLawPtr(CL_LE_ID);
    auto CL_SCSB_Ptr = S.ConstitutiveLawGetConstitutiveLawPtr(CL_SCSB_ID);
    auto CL_MT_Ptr = S.ConstitutiveLawGetConstitutiveLawPtr(CL_MT_ID);
    auto CL_AIE_Ptr = static_cast<NuTo::AdditiveInputExplicit*>(S.ConstitutiveLawGetConstitutiveLawPtr(CL_AIE_ID));
    auto CL_AO_Ptr = static_cast<NuTo::AdditiveOutput*>(S.ConstitutiveLawGetConstitutiveLawPtr(CL_AO_ID));

    TimeControl tCtrl;
    tCtrl.t_final = 365.0 * 24.0 * 60.0 * 60.0;
    tCtrl.delta_t = tCtrl.t_final;
    tCtrl.t_write = tCtrl.t_final;


    MoistureTransportControl MTCtrl(S, *CL_MT_Ptr);
    MTCtrl.InitialRelativeHumidity = 1.0;
    MTCtrl.BoundaryEnvironmentalRH = 0.40;
    MTCtrl.MassExchangeRate = 1.;
    MTCtrl.DiffusionCoefficientRH = 1.e-5;
    MTCtrl.DiffusionCoefficientWV = 1.e-1;
    MTCtrl.BoundaryDiffusionCoefficientRH = 1.0e-5;
    MTCtrl.BoundaryDiffusionCoefficientWV = 1.0e-1;
    MTCtrl.SetParametersConstitutiveLaw();

    MechanicsControl MeCtrl(S, *CL_LE_Ptr);
    MeCtrl.SetParametersConstitutiveLaw();

    CL_SCSB_Ptr->SetParameterDouble(NuTo::Constitutive::eConstitutiveParameter::MACROSCOPIC_BULK_MODULUS,
                                    TEST_MACROSCOPICBULKMODULUS);
    CL_SCSB_Ptr->SetParameterDouble(NuTo::Constitutive::eConstitutiveParameter::SOLID_PHASE_BULK_MODULUS,
                                    TEST_SOLIDPHASEBULKMODULUS);
    CL_SCSB_Ptr->SetParameterDouble(NuTo::Constitutive::eConstitutiveParameter::TEMPERATURE, TEST_TEMPERATURE);


    CL_AIE_Ptr->AddConstitutiveLaw(*CL_LE_Ptr);
    CL_AIE_Ptr->AddConstitutiveLaw(*CL_SCSB_Ptr, NuTo::Constitutive::eInput::ENGINEERING_STRAIN);

    CL_AO_Ptr->AddConstitutiveLaw(*CL_MT_Ptr);
    CL_AO_Ptr->AddConstitutiveLaw(*CL_AIE_Ptr);

    SetupStructure(S, testName);
    auto SEC = SetupSection<TDim>(S);

    auto meshInfo = NuTo::MeshGenerator::Grid(S, rL, rN);

    for (auto& it : rDofIPTMap)
        S.InterpolationTypeAdd(meshInfo.second, it.first, it.second);

    S.ElementGroupSetSection(meshInfo.first, SEC);
    S.ElementGroupSetConstitutiveLaw(meshInfo.first, CL_AO_ID);

    SetupIntegrationType<TDim>(S, meshInfo.first);

    S.ElementTotalConvertToInterpolationType(); // old used values 1.0e-12,0.001
    MTCtrl.ApplyInitialNodalValues();


    auto LambdaGetBoundaryNodes = [rL](NuTo::NodeBase* rNodePtr) -> bool {
        double Tol = 1.e-6;
        if (rNodePtr->GetNum(NuTo::Node::eDof::COORDINATES) > 0)
        {
            double x = rNodePtr->Get(NuTo::Node::eDof::COORDINATES)[0];
            if ((x >= 0.0 - Tol && x <= 0.0 + Tol) || (x >= rL[0] - Tol && x <= rL[0] + Tol))
            {
                return true;
            }

            if (TDim > 1)
            {
                double y = rNodePtr->Get(NuTo::Node::eDof::COORDINATES)[1];
                if ((y >= 0.0 - Tol && y <= 0.0 + Tol) || (y >= rL[1] - Tol && y <= rL[1] + Tol))
                {
                    return true;
                }
            }
            if (TDim > 2)
            {
                double z = rNodePtr->Get(NuTo::Node::eDof::COORDINATES)[2];
                if ((z >= 0.0 - Tol && z <= 0.0 + Tol) || (z >= rL[2] - Tol && z <= rL[2] + Tol))
                {
                    return true;
                }
            }
        }
        return false;
    }; // GetBoundaryNodesLambda


    auto LambdaTimeDepBoundaryRH = [&MTCtrl, &tCtrl](double rTime) -> double {
        if (rTime == 0.0)
        {
            return MTCtrl.InitialRelativeHumidity;
        }
        else
        {
            if (rTime < tCtrl.BC_TransitionTime)
            {
                return MTCtrl.InitialRelativeHumidity -
                       sin(rTime / tCtrl.BC_TransitionTime * 3.14 / 2.0) *
                               (MTCtrl.InitialRelativeHumidity - MTCtrl.BoundaryEnvironmentalRH);
            }
            {
                return MTCtrl.BoundaryEnvironmentalRH;
            }
        }
    }; // TimeDepBoundaryRHLambda


    SetupConstrainedNodeBoundaryElements<TDim>(S, LambdaGetBoundaryNodes, LambdaTimeDepBoundaryRH);


    auto lambdaGetNodeLeftBottomFront = [rL](NuTo::NodeBase* rNodePtr) -> bool {
        if (rNodePtr->GetNum(NuTo::Node::eDof::DISPLACEMENTS) == 0)
            return false;
        double Tol = 1.e-6;
        if (rNodePtr->GetNum(NuTo::Node::eDof::COORDINATES) > 0)
        {
            double x = 0.0, y = 0.0, z = 0.0;
            x = rNodePtr->Get(NuTo::Node::eDof::COORDINATES)[0];
            if (TDim > 1)
                y = rNodePtr->Get(NuTo::Node::eDof::COORDINATES)[1];
            if (TDim > 2)
                z = rNodePtr->Get(NuTo::Node::eDof::COORDINATES)[2];

            if (x >= 0.0 - Tol && x <= 0.0 + Tol && y >= 0.0 - Tol && y <= 0.0 + Tol && z >= 0.0 - Tol &&
                z <= 0.0 + Tol)
            {
                return true;
            }
        }
        return false;
    }; // lambdaGetNodeLeftBottom

    auto lambdaGetNodeLeftTopFront = [rL](NuTo::NodeBase* rNodePtr) -> bool {
        if (rNodePtr->GetNum(NuTo::Node::eDof::DISPLACEMENTS) == 0)
            return false;
        double Tol = 1.e-6;
        if (rNodePtr->GetNum(NuTo::Node::eDof::COORDINATES) > 0)
        {
            double x = 0.0, y = 0.0, z = 0.0;
            x = rNodePtr->Get(NuTo::Node::eDof::COORDINATES)[0];
            if (TDim > 1)
                y = rNodePtr->Get(NuTo::Node::eDof::COORDINATES)[1];
            if (TDim > 2)
                z = rNodePtr->Get(NuTo::Node::eDof::COORDINATES)[2];

            if (x >= 0.0 - Tol && x <= 0.0 + Tol && y >= 0.0 - Tol && y <= 0.0 + Tol && z >= rL[2] - Tol &&
                z <= rL[2] + Tol)
            {
                return true;
            }
        }
        return false;
    }; // lambdaGetNodeLeftBottom

    auto lambdaGetNodeLeftBottomBack = [rL](NuTo::NodeBase* rNodePtr) -> bool {
        if (rNodePtr->GetNum(NuTo::Node::eDof::DISPLACEMENTS) == 0)
            return false;
        double Tol = 1.e-6;
        if (rNodePtr->GetNum(NuTo::Node::eDof::COORDINATES) > 0)
        {
            double x = 0.0, y = 0.0, z = 0.0;
            x = rNodePtr->Get(NuTo::Node::eDof::COORDINATES)[0];
            if (TDim > 1)
                y = rNodePtr->Get(NuTo::Node::eDof::COORDINATES)[1];
            if (TDim > 2)
                z = rNodePtr->Get(NuTo::Node::eDof::COORDINATES)[2];

            if (x >= 0.0 - Tol && x <= 0.0 + Tol && y >= rL[1] - Tol && y <= rL[1] + Tol && z >= 0.0 - Tol &&
                z <= 0.0 + Tol)
            {
                return true;
            }
        }
        return false;
    }; // lambdaGetNodeLeftBottom


    MeCtrl.AddConstraint(lambdaGetNodeLeftBottomFront, NuTo::eDirection::X);
    if (TDim > 1)
    {
        MeCtrl.AddConstraint(lambdaGetNodeLeftBottomFront, NuTo::eDirection::Y);
        MeCtrl.AddConstraint(lambdaGetNodeLeftBottomBack, NuTo::eDirection::X);
    }

    if (TDim > 2)
    {
        MeCtrl.AddConstraint(lambdaGetNodeLeftBottomFront, NuTo::eDirection::Z);
        MeCtrl.AddConstraint(lambdaGetNodeLeftBottomBack, NuTo::eDirection::Z);
        MeCtrl.AddConstraint(lambdaGetNodeLeftTopFront, NuTo::eDirection::X);
        MeCtrl.AddConstraint(lambdaGetNodeLeftTopFront, NuTo::eDirection::Y);
    }


    MTCtrl.SetupStaticData();
    S.NodeBuildGlobalDofs();

    SetupMultiProcessor(S);

    SetupVisualize(S, true);

    SetupTimeIntegration(TI, tCtrl, resultDir, rStaggered);
    NuTo::Timer timer("shrinkagetest");
    TI.Solve(tCtrl.t_final);

    CheckMoistureTransportResults<TDim>(S, rN, rL);
    /*CheckMechanicsResults<TDim>(S)*/;
}

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//  Main
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

int main()
{
    std::map<NuTo::Node::eDof, NuTo::Interpolation::eTypeOrder> dofIPTMap;
    dofIPTMap[NuTo::Node::eDof::DISPLACEMENTS] = NuTo::Interpolation::eTypeOrder::EQUIDISTANT1;
    dofIPTMap[NuTo::Node::eDof::RELATIVEHUMIDITY] = NuTo::Interpolation::eTypeOrder::EQUIDISTANT1;
    dofIPTMap[NuTo::Node::eDof::WATERVOLUMEFRACTION] = NuTo::Interpolation::eTypeOrder::EQUIDISTANT1;


    // STRESS based

    ShrinkageTestStressBased<1>({3}, {0.01}, dofIPTMap);

    ShrinkageTestStressBased<2>({3, 3}, {0.01, 0.01}, dofIPTMap);

    ShrinkageTestStressBased<3>({3, 3, 3}, {0.01, 0.01, 0.01}, dofIPTMap);


    ShrinkageTestStressBased<1>({3}, {0.01}, dofIPTMap, true);

    ShrinkageTestStressBased<2>({3, 3}, {0.01, 0.01}, dofIPTMap, true);

    ShrinkageTestStressBased<3>({3, 3, 3}, {0.01, 0.01, 0.01}, dofIPTMap, true);

    // STRAIN based

    ShrinkageTestStrainBased<1>({3}, {0.01}, dofIPTMap);

    ShrinkageTestStrainBased<2>({3, 3}, {0.01, 0.01}, dofIPTMap);

    ShrinkageTestStrainBased<3>({3, 3, 3}, {0.01, 0.01, 0.01}, dofIPTMap);


    ShrinkageTestStrainBased<1>({3}, {0.01}, dofIPTMap, true);

    ShrinkageTestStrainBased<2>({3, 3}, {0.01, 0.01}, dofIPTMap, true);

    ShrinkageTestStrainBased<3>({3, 3, 3}, {0.01, 0.01, 0.01}, dofIPTMap, true);
}
