#include "nuto/base/Timer.h"
#include "nuto/mechanics/structures/unstructured/Structure.h"
#include "nuto/mechanics/timeIntegration/NewmarkDirect.h"
#include "nuto/mechanics/tools/MeshGenerator.h"
#include "nuto/physics/PhysicalConstantsSI.h"
#include "nuto/physics/PhysicalEquationsSI.h"
#include <array>


#include "nuto/mechanics/constitutive/laws/MoistureTransport.h"


//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//  Setup Preprocessor
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


/*---------------------------------------------*\
|*                  DEFINES                    *|
\*---------------------------------------------*/


// --- Time integration scheme
// ---------------------------
#define RES_TOLERANCE_MECHANICS 1e-4
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
    double          delta_t                         = 1.0/1.0  *     1.0 * 24.0 * 60.0 * 60.0;
    double          t_write                         = 1.0/1.0  *     1.0 * 24.0 * 60.0 * 60.0;
    double          t_final                         = 20.0/1.0 *     1.0 * 24.0 * 60.0 * 60.0;
    double          BC_TransitionTime               =                      24.0 * 60.0 * 60.0;
};

struct MechanicsControl
{
    //! @brief constructor
    MechanicsControl( NuTo::Structure& rS,
                      NuTo::ConstitutiveBase& rCL)
        : mS(rS),
          mCL(rCL)
    {
        switch(mCL.GetType())
        {
        case NuTo::Constitutive::eConstitutiveType::LINEAR_ELASTIC_ENGINEERING_STRESS:
            break;
        default:
            throw NuTo::Exception(__PRETTY_FUNCTION__,"The referenced constitutive law is not a mechanics model");
        }
    }


    void SetParametersConstitutiveLaw()
    {
        // set variables

        mCL.SetParameterDouble  (NuTo::Constitutive::eConstitutiveParameter::DENSITY,           mDensity);
        mCL.SetParameterDouble  (NuTo::Constitutive::eConstitutiveParameter::POISSONS_RATIO,    mPoissonRatio);
        mCL.SetParameterDouble  (NuTo::Constitutive::eConstitutiveParameter::YOUNGS_MODULUS,    mYoungsModulus);
    }

    template<int TDim>
    int AddConstraint(NuTo::NewmarkDirect& rTI,
                       std::function<bool(NuTo::NodeBase*)> rGetNodeFunction,
                       unsigned int rDirection,
                       double rValue = 0.0)
    {
        assert(rValue <= TDim && "Direction isn't part of current dimension");
        int GRPNodesConstraint = mS.GroupCreate("Nodes");
        mS.GroupAddNodeFunction(GRPNodesConstraint,rGetNodeFunction);

        NuTo::FullMatrix<double,Eigen::Dynamic,Eigen::Dynamic> direction(TDim,1);
        direction.SetValue(rDirection, 0, 1.0);

        return mS.ConstraintLinearSetDisplacementNodeGroup(GRPNodesConstraint, direction,rValue);


    }

    template<int TDim>
    int AddTimeDependentConstraint( NuTo::NewmarkDirect& rTI,
                                    std::function<bool(NuTo::NodeBase*)> rGetNodeFunction,
                                    std::function<double(double)> rDisplacementFunction,
                                    unsigned int rDirection,
                                    double rValue = 0.0)
    {
        int constraintID = AddConstraint<TDim>(rTI,
                                               rGetNodeFunction,
                                               rDirection,
                                               rValue);
        rTI.AddTimeDependentConstraintFunction(constraintID,
                                               rDisplacementFunction);

        return constraintID;
    }


    NuTo::Structure&        mS;
    NuTo::ConstitutiveBase& mCL;

    double mYoungsModulus   =   TEST_YOUNGSMODULUS;
    double mPoissonRatio    =   0.0;
    double mDensity         =   1.0;
};


struct MoistureTransportControl
{

    // ctor
    MoistureTransportControl(NuTo::Structure& rS,
                             NuTo::ConstitutiveBase& rMT)
        : mS(rS),
          mMT(rMT)
    {
        if(mMT.GetType()!=NuTo::Constitutive::eConstitutiveType::MOISTURE_TRANSPORT)
            throw NuTo::Exception(__PRETTY_FUNCTION__,"The referenced constitutive law is not a moisture transport model");

        // Values fitted from figure in Johannessons paper
        AdsorptionCoeffs(0) =  0.0;
        AdsorptionCoeffs(1) =  0.19692057340725558;
        AdsorptionCoeffs(2) = -0.28253538941816925;
        AdsorptionCoeffs(3) =  0.22661481601091368;

        DesorptionCoeffs(0) =  0.0;
        DesorptionCoeffs(1) =  0.26719233184420238;
        DesorptionCoeffs(2) = -0.41030868184510738;
        DesorptionCoeffs(3) =  0.32511635000090505;
    }

public:
    void SetParametersConstitutiveLaw()
    {
        // set variables
        mMT.SetParameterBool    (NuTo::Constitutive::eConstitutiveParameter::ENABLE_MODIFIED_TANGENTIAL_STIFFNESS,      EnableModiefiedTangentialStiffness);
        mMT.SetParameterBool    (NuTo::Constitutive::eConstitutiveParameter::ENABLE_SORPTION_HYSTERESIS,                EnableSorptionHysteresis);

        mMT.SetParameterDouble  (NuTo::Constitutive::eConstitutiveParameter::BOUNDARY_DIFFUSION_COEFFICIENT_RH,         BoundaryDiffusionCoefficientRH);
        mMT.SetParameterDouble  (NuTo::Constitutive::eConstitutiveParameter::BOUNDARY_DIFFUSION_COEFFICIENT_WV,         BoundaryDiffusionCoefficientWV);
        mMT.SetParameterDouble  (NuTo::Constitutive::eConstitutiveParameter::DENSITY_WATER,                             DensityWater);
        mMT.SetParameterDouble  (NuTo::Constitutive::eConstitutiveParameter::DIFFUSION_COEFFICIENT_RH,                  DiffusionCoefficientRH);
        mMT.SetParameterDouble  (NuTo::Constitutive::eConstitutiveParameter::DIFFUSION_COEFFICIENT_WV,                  DiffusionCoefficientWV);
        mMT.SetParameterDouble  (NuTo::Constitutive::eConstitutiveParameter::DIFFUSION_EXPONENT_RH,                     DiffusionExponentRH);
        mMT.SetParameterDouble  (NuTo::Constitutive::eConstitutiveParameter::DIFFUSION_EXPONENT_WV,                     DiffusionExponentWV);
        mMT.SetParameterDouble  (NuTo::Constitutive::eConstitutiveParameter::GRADIENT_CORRECTION_ADSORPTION_DESORPTION, GradientCorrectionAdsorpionDesorption);
        mMT.SetParameterDouble  (NuTo::Constitutive::eConstitutiveParameter::GRADIENT_CORRECTION_DESORPTION_ADSORPTION, GradientCorrectionDesorptionAdsorpion);
        mMT.SetParameterDouble  (NuTo::Constitutive::eConstitutiveParameter::MASS_EXCHANGE_RATE,                        MassExchangeRate);
        mMT.SetParameterDouble  (NuTo::Constitutive::eConstitutiveParameter::PORE_VOLUME_FRACTION,                      PoreVolumeFraction);
        mMT.SetParameterDouble  (NuTo::Constitutive::eConstitutiveParameter::DENSITY_SATURATED_WATER_VAPOR,             DensitySaturatedWaterVapor);

        mMT.SetParameterFullVectorDouble    (NuTo::Constitutive::eConstitutiveParameter::POLYNOMIAL_COEFFICIENTS_ADSORPTION,    AdsorptionCoeffs);
        mMT.SetParameterFullVectorDouble    (NuTo::Constitutive::eConstitutiveParameter::POLYNOMIAL_COEFFICIENTS_DESORPTION,    DesorptionCoeffs);


        // Calculate equilibrium water volume fraction
        InitialWaterVolumeFraction   = mMT.GetEquilibriumWaterVolumeFraction(InitialRelativeHumidity,
                                                                             DesorptionCoeffs);
    }



public:
    void ApplyInitialNodalValues()
    {
        unsigned int NNodes = mS.GetNumNodes();

        for (unsigned int i=0; i<NNodes; i++)
        {


            if(mS.NodeGetNodePtr(i)->GetNumRelativeHumidity() != 0)
            {
                mS.NodeGetNodePtr(i)->SetRelativeHumidity(0,InitialRelativeHumidity) ;
            }
            if(mS.NodeGetNodePtr(i)->GetNumWaterVolumeFraction() != 0)
            {
                mS.NodeGetNodePtr(i)->SetWaterVolumeFraction(0,InitialWaterVolumeFraction);
            }
        }
    }

    void SetupStaticData()
    {
        for (int i=0; i<mS.GetNumElements(); i++)
        {
            for (int theIP=0; theIP< mS.ElementGetElementPtr(i)->GetNumIntegrationPoints(); theIP++)
            {
                NuTo::ConstitutiveStaticDataMoistureTransport *StaticData = mS.ElementGetElementPtr(i)->GetStaticData(theIP)->AsMoistureTransport();
                StaticData->SetLastSorptionCoeff(mMT.GetParameterFullVectorDouble(NuTo::Constitutive::eConstitutiveParameter::POLYNOMIAL_COEFFICIENTS_DESORPTION));
                StaticData->SetCurrentSorptionCoeff(mMT.GetParameterFullVectorDouble(NuTo::Constitutive::eConstitutiveParameter::POLYNOMIAL_COEFFICIENTS_DESORPTION));
                StaticData->SetLastRelHumValue(InitialRelativeHumidity);
                StaticData->SetSorptionHistoryDesorption(SorptionHistoryDesorption);
            }
        }
    }




    // general
    bool            EnableModiefiedTangentialStiffness  =   false;

    double          InitialRelativeHumidity                 =   1.0;
    double          InitialWaterVolumeFraction              =   0.;         //! --> Calculated from relative humidity
    double          MassExchangeRate                        =   3.42e-7;
    double          PoreVolumeFraction                      =   0.25;
    double          DiffusionCoefficientRH                  =   3.9e-12;
    double          DiffusionExponentRH                     =   1.0;
    double          DensitySaturatedWaterVapor              =   0.0173;
    double          DensityWater                            =   999.97;
    double          DiffusionCoefficientWV                  =   1.17e-7;
    double          DiffusionExponentWV                     =   2.0;

    // boundary condition
    double          BoundaryEnvironmentalRH                 =   0.45;
    double          BoundaryDiffusionCoefficientRH          =   1.0e-10 * 1000;
    double          BoundaryDiffusionCoefficientWV          =   1.0e-7 * 1000;

    // sorption isotherms
    bool            EnableSorptionHysteresis                =   false;
    bool            SorptionHistoryDesorption               =   true;
    double          GradientCorrectionDesorptionAdsorpion   =   0.26;
    double          GradientCorrectionAdsorpionDesorption   =   0.56;
    NuTo::FullVector<double,4> AdsorptionCoeffs             = Eigen::MatrixXd::Constant(4, 1, 0.0);
    NuTo::FullVector<double,4> DesorptionCoeffs             = Eigen::MatrixXd::Constant(4, 1, 0.0);


    //references
    NuTo::Structure&            mS;
    NuTo::ConstitutiveBase&     mMT;
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
                                          NuTo::NewmarkDirect& rTI,
                                          std::function<double(double)> rBoundaryConstraintFunction)
{
    int nGrpBE = rS.GroupCreate("NODES");
    rS.GroupAddNodeFunction(nGrpBE,rFunctionGetBoundaryNode);

    int eGrpBE = rS.GroupCreate("ELEMENTS");
    rS.GroupAddElementsFromNodes(eGrpBE, nGrpBE, false);

    std::set<NuTo::Node::eDof> controlNodeDofs;
    controlNodeDofs.insert(NuTo::Node::RELATIVEHUMIDITY);

    int boundaryControlNodeID = rS.NodeCreateDOFs(controlNodeDofs);
    NuTo::NodeBase* controlNodePtr = rS.NodeGetNodePtr(boundaryControlNodeID);
    int groupBoundaryElements = rS.BoundaryElementsCreate(eGrpBE,nGrpBE,controlNodePtr);
    int controlNodeConstraint = rS.ConstraintLinearSetRelativeHumidityNode(controlNodePtr,1.0);

    // Set Integration type - default not sufficient
    NuTo::FullVector<int,Eigen::Dynamic> boundaryElementIDs;
    rS.ElementGroupGetMembers(groupBoundaryElements, boundaryElementIDs);


    for(unsigned int i=0; i<boundaryElementIDs.rows();++i)
    {
        NuTo::ElementBase* elementPtr =  rS.ElementGetElementPtr(boundaryElementIDs[i]);
        switch(TDim)
        {
        case 1:
            elementPtr->SetIntegrationType(rS.GetPtrIntegrationType(NuTo::IntegrationType::IntegrationType0DBoundary), elementPtr->GetIpDataType(0));
            break;
        case 2:
            elementPtr->SetIntegrationType(rS.GetPtrIntegrationType(NuTo::IntegrationType::IntegrationType1D2NGauss2Ip), elementPtr->GetIpDataType(0));
            break;
        case 3:
            elementPtr->SetIntegrationType(rS.GetPtrIntegrationType(NuTo::IntegrationType::IntegrationType2D4NGauss4Ip), elementPtr->GetIpDataType(0));
            break;
        default:
            throw NuTo::Exception(__PRETTY_FUNCTION__,"Invalid dimension");
        }


    }

    rTI.AddTimeDependentConstraintFunction(controlNodeConstraint, rBoundaryConstraintFunction);
}


/*---------------------------------------------*\
|*              integration type               *|
\*---------------------------------------------*/

template <int TDim>
void SetupIntegrationType(NuTo::Structure& rS, int rIPT)
{
    switch(TDim)
    {
    case 1:
        rS.InterpolationTypeSetIntegrationType(rIPT,NuTo::IntegrationType::IntegrationType1D2NGauss2Ip,NuTo::IpData::STATICDATA);
        break;
    case 2:
        rS.InterpolationTypeSetIntegrationType(rIPT,NuTo::IntegrationType::IntegrationType2D4NGauss4Ip,NuTo::IpData::STATICDATA);
        break;
    case 3:
        rS.InterpolationTypeSetIntegrationType(rIPT,NuTo::IntegrationType::IntegrationType3D8NGauss2x2x2Ip,NuTo::IpData::STATICDATA);
        break;
    default:
        throw NuTo::Exception(__PRETTY_FUNCTION__,"Invalid dimension");
    }

}


/*---------------------------------------------*\
|*             interpolation type              *|
\*---------------------------------------------*/

template <int TDim>
inline int SetupInterpolationType(NuTo::Structure& rS,
                                  std::map<NuTo::Node::eDof,NuTo::Interpolation::eTypeOrder> rDofIPTMap,
                                  std::string rShape = "")
{
    if(rShape.empty())
    {
        switch(TDim)
        {
        case 1:
            rShape = "TRUSS1D";
            break;

        case 2:
            rShape = "QUAD2D";
            break;

        case 3:
            rShape = "BRICK3D";
            break;

        default:
            throw NuTo::Exception(__PRETTY_FUNCTION__,"Invalid dimension");
        }
    }


    int IPT = rS.InterpolationTypeCreate(rShape);
    for(auto itIPT : rDofIPTMap)
    {
        rS.InterpolationTypeAdd(IPT, itIPT.first, itIPT.second);
    }
    return IPT;
}

/*---------------------------------------------*\
|*                  mesh setup                 *|
\*---------------------------------------------*/

template<int TDim>
void SetupMesh(NuTo::Structure &rS,
               int rSEC,
               int rConsLaw,
               int rIPT,
               std::array<int, TDim> rN,
               std::array<double,TDim> rL)
{
    throw NuTo::Exception(__PRETTY_FUNCTION__,"Invalid dimension");
}



template<>
void SetupMesh<1>(NuTo::Structure &rS,
                  int rSEC,
                  int rConsLaw,
                  int rIPT,
                  std::array<int, 1> rN,
                  std::array<double,1> rL)
{
NuTo::MeshGenerator::MeshLineSegment(rS,
                                     rSEC,
                                     rConsLaw,
                                     rIPT,
                                     rN,
                                     rL);
}

template<>
void SetupMesh<2>(NuTo::Structure &rS,
                  int rSEC,
                  int rConsLaw,
                  int rIPT,
                  std::array<int, 2> rN,
                  std::array<double,2> rL)
{
NuTo::MeshGenerator::MeshRectangularPlane(rS,
                                          rSEC,
                                          rConsLaw,
                                          rIPT,
                                          rN,
                                          rL);
}


template<>
void SetupMesh<3>(NuTo::Structure &rS,
                  int rSEC,
                  int rConsLaw,
                  int rIPT,
                  std::array<int, 3> rN,
                  std::array<double,3> rL)
{
NuTo::MeshGenerator::MeshCuboid(rS,
                                rSEC,
                                rConsLaw,
                                rIPT,
                                rN,
                                rL);

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
inline int SetupSection(NuTo::Structure& rS, double rAreaThickness = 1.0)
{
    switch (TDim)
    {
    case 1:
    {
        int Sec = rS.SectionCreate("TRUSS");
        rS.SectionSetArea(Sec,rAreaThickness);
        rS.ElementTotalSetSection(Sec);
        return Sec;
    }

    case 2:
    {
        int Sec = rS.SectionCreate("PLANE_STRESS");
        rS.SectionSetThickness(Sec,rAreaThickness);
        rS.ElementTotalSetSection(Sec);
        return Sec;
    }

    case 3:
    {
        int Sec = rS.SectionCreate("VOLUME");
        rS.ElementTotalSetSection(Sec);
        return Sec;
    }

    default:
        throw NuTo::Exception(__PRETTY_FUNCTION__,"Invalid dimension");
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
    Log.OpenFile(rTestName+".log");
}


/*---------------------------------------------*\
|*              time integration               *|
\*---------------------------------------------*/

inline void SetupTimeIntegration(NuTo::NewmarkDirect& rTI,
                                 const TimeControl& rTC,
                                 const std::string& rResultDir,
                                 bool rStaggered)
{
    rTI.SetPerformLineSearch(false);
    rTI.SetVerboseLevel(0);
    rTI.SetToleranceResidual(NuTo::Node::eDof::DISPLACEMENTS,RES_TOLERANCE_MECHANICS);
    rTI.SetToleranceResidual(NuTo::Node::eDof::RELATIVEHUMIDITY,RES_TOLERANCE_MOISTURE_TRANSPORT);
    rTI.SetToleranceResidual(NuTo::Node::eDof::WATERVOLUMEFRACTION,RES_TOLERANCE_MOISTURE_TRANSPORT);
    rTI.SetMaxNumIterations(MAX_ITERATION);

    rTI.SetTimeStep(rTC.delta_t);
    rTI.SetMinTimeStepPlot(rTC.t_write);

    rTI.SetResultDirectory(rResultDir,true);

    if(rStaggered)
    {
        rTI.AddCalculationStep({NuTo::Node::eDof::RELATIVEHUMIDITY,NuTo::Node::eDof::WATERVOLUMEFRACTION});
        rTI.AddCalculationStep({NuTo::Node::eDof::DISPLACEMENTS});
    }
}

/*---------------------------------------------*\
|*                 visualize                   *|
\*---------------------------------------------*/

inline void SetupVisualize(NuTo::Structure& rS)
{
#ifdef ENABLE_VISUALIZE
        int visGrp = rS.GroupCreate(NuTo::Groups::eGroupId::Elements);
        rS.GroupAddElementsTotal(visGrp);
        rS.AddVisualizationComponent(visGrp, NuTo::VisualizeBase::DISPLACEMENTS);
        rS.AddVisualizationComponent(visGrp, NuTo::VisualizeBase::RELATIVE_HUMIDITY);
        rS.AddVisualizationComponent(visGrp, NuTo::VisualizeBase::WATER_VOLUME_FRACTION);
        rS.AddVisualizationComponent(visGrp, NuTo::VisualizeBase::ENGINEERING_STRAIN);
        rS.AddVisualizationComponent(visGrp, NuTo::VisualizeBase::PRINCIPAL_ENGINEERING_STRESS);
#endif // ENABLE_VISUALIZE
}








//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//  Check Results
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

template<int TDim>
void CheckMechanicsResults(NuTo::Structure& rS)
{
    constexpr const double temperature     =   TEST_TEMPERATURE;
    constexpr const double capStressFactor =   NuTo::SI::DensityLiquidWater(temperature) * NuTo::SI::IdealGasConstant
                                             * temperature / NuTo::SI::MolarMassWater;

    constexpr const double capStress       = 0.062035 * log(0.4) * capStressFactor;

    const NodeMap& nodePtrMap = rS.NodeGetNodeMap();
    BOOST_FOREACH(NodeMap::const_iterator::value_type it, nodePtrMap)
    {
        const NuTo::NodeBase* nodePtr = it.second;
        if(nodePtr->GetNumDisplacements()<1)
        {
            continue;   // Nodes without Displacements cant be checked
        }
        for(int i=0; i<TDim; ++i)
        {
            double coord = nodePtr->GetCoordinate(i);
            double refDisp = capStress / TEST_YOUNGSMODULUS * coord;
            double disp  = nodePtr->GetDisplacement(i);
            double diff   = std::abs(refDisp) - std::abs(disp);
            constexpr const double tolerance = 1e-10;
            if((diff >tolerance || diff < -tolerance))// && coordX > 0)
            {
                throw NuTo::Exception(__PRETTY_FUNCTION__,"One ore more calculated Displacements are not correct");
            }
        }
    }
        std::cout << "Displacements correct!" << std::endl;
}

template<int TDim>
void CheckMoistureTransportResults(NuTo::Structure& rS,
                                   std::array<int,TDim> rN,
                                   std::array<double,TDim> rL)
{


    constexpr const double tolerance = 0.00001; // Tolerance because not all necessary value (sorption curve) are given in the paper and must be approximated
    unsigned int numMismatchingValues = 0;


    const NodeMap& nodePtrMap = rS.NodeGetNodeMap();
    BOOST_FOREACH(NodeMap::const_iterator::value_type it, nodePtrMap)
    {
        const NuTo::NodeBase* nodePtr = it.second;


        if(nodePtr->GetNumWaterVolumeFraction()>0)
        {
            double nodalWVF = nodePtr->GetWaterVolumeFraction();
            double eqWVF = 0.062035;
            if(nodalWVF<eqWVF-tolerance || nodalWVF>eqWVF+tolerance)
            {
                ++numMismatchingValues;
            }
        }
        if(nodePtr->GetNumWaterVolumeFraction()>0)
        {
            double nodalRH = nodePtr->GetRelativeHumidity();
            double eqRH = 0.4;
            if(nodalRH<eqRH-tolerance || nodalRH>eqRH+tolerance)
            {
                ++numMismatchingValues;
            }
        }


    }
    if(numMismatchingValues>0)
    {
        throw NuTo::Exception(__PRETTY_FUNCTION__,"One ore more calculated relative humidity/water volue fraction values exceeds the tolerance when compared to reference values");
    }
    std::cout << "Water volume fraction and relative humidity correct!" << std::endl;
}



//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//  Tests
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



//! @brief performs a simulation in the desired dimension
//! @param rN: array with number of elements in each direction
//! @param rL: array with length of elements in each direction
template<int TDim>
void ShrinkageTest(std::array<int,TDim> rN,
                        std::array<double,TDim> rL,
                        std::map<NuTo::Node::eDof,NuTo::Interpolation::eTypeOrder> rDofIPTMap,
                        bool rStaggered = false)
{
//    assert((rN[0] == 16)   && "Only 16 elements in x-direction allowed for this test --- needed for direct comparison to values given in Johannesson 2010");
//    assert((rL[0] == 0.16) && "The length in flow direction (x) must be 0.16m for direct comparison with paper values");


    std::string testName = std::string("ShrinkageTest") + std::to_string(TDim) +"D";
    if(rStaggered)
        testName += "_staggered";
    std::string resultDir = std::string("./MultipleConstitutiveLaws_") + testName;

    std::cout << std::endl << "--------------------------------------------------------------------------"
              << std::endl << "Start test: "<< testName
              << std::endl << "--------------------------------------------------------------------------" << std::endl;

    // Allocate neccessary stuff
    NuTo::Structure S(TDim);
    NuTo::NewmarkDirect TI(&S);
    int CL_LE_ID   = S.ConstitutiveLawCreate(NuTo::Constitutive::eConstitutiveType::LINEAR_ELASTIC_ENGINEERING_STRESS);
    int CL_SCSB_ID = S.ConstitutiveLawCreate(NuTo::Constitutive::eConstitutiveType::SHRINKAGE_CAPILLARY_STRESS_BASED);
    int CL_MT_ID   = S.ConstitutiveLawCreate(NuTo::Constitutive::eConstitutiveType::MOISTURE_TRANSPORT);
    int CL_CLAL_ID = S.ConstitutiveLawCreate(NuTo::Constitutive::eConstitutiveType::CONSTITUTIVE_LAWS_ADDITIVE_OUTPUT);



    NuTo::ConstitutiveBase* CL_LE_Ptr   = S.ConstitutiveLawGetConstitutiveLawPtr(CL_LE_ID);
    NuTo::ConstitutiveBase* CL_SCSB_Ptr = S.ConstitutiveLawGetConstitutiveLawPtr(CL_SCSB_ID);
    NuTo::ConstitutiveBase* CL_MT_Ptr   = S.ConstitutiveLawGetConstitutiveLawPtr(CL_MT_ID);
    NuTo::ConstitutiveBase* CL_CLAL_Ptr = S.ConstitutiveLawGetConstitutiveLawPtr(CL_CLAL_ID);

    TimeControl tCtrl;
    tCtrl.t_final = 365.0 * 24.0 * 60.0 * 60.0;
    tCtrl.delta_t = tCtrl.t_final;
    tCtrl.t_write = tCtrl.t_final;


    MoistureTransportControl MTCtrl(S,*CL_MT_Ptr);
    MTCtrl.InitialRelativeHumidity = 0.95;
    MTCtrl.BoundaryEnvironmentalRH = 0.40;
    MTCtrl.MassExchangeRate = 1.;
    MTCtrl.DiffusionCoefficientRH = 1e-5;
    MTCtrl.DiffusionCoefficientWV = 1.;
    MTCtrl.BoundaryDiffusionCoefficientRH          =   1.0e-4;
    MTCtrl.BoundaryDiffusionCoefficientWV          =   1.0;
    MTCtrl.SetParametersConstitutiveLaw();

    MechanicsControl MeCtrl(S,*CL_LE_Ptr);
    MeCtrl.SetParametersConstitutiveLaw();

    CL_SCSB_Ptr->SetParameterDouble(NuTo::Constitutive::eConstitutiveParameter::TEMPERATURE, TEST_TEMPERATURE);


    CL_CLAL_Ptr->AddConstitutiveLaw(CL_LE_Ptr);
    CL_CLAL_Ptr->AddConstitutiveLaw(CL_SCSB_Ptr);
    CL_CLAL_Ptr->AddConstitutiveLaw(CL_MT_Ptr);



    SetupStructure(S,testName);
    int SEC = SetupSection<TDim>(S);
    int IPT = SetupInterpolationType<TDim>(S,rDofIPTMap);

    SetupMesh<TDim>(S,
                    SEC,
                    CL_CLAL_ID,
                    IPT,
                    rN,
                    rL);

    SetupIntegrationType<TDim>(S,IPT);

    S.ElementTotalConvertToInterpolationType(); //old used values 1.0e-12,0.001
    MTCtrl.ApplyInitialNodalValues();


    auto LambdaGetBoundaryNodes = [rL](NuTo::NodeBase* rNodePtr) -> bool
                                {
                                    double Tol = 1.e-6;
                                    if (rNodePtr->GetNumCoordinates()>0)
                                    {
                                        double x = rNodePtr->GetCoordinate(0);
                                        if ((x >= 0.0   - Tol   && x <= 0.0   + Tol) ||
                                            (x >= rL[0] - Tol   && x <= rL[0] + Tol))
                                        {
                                            return true;
                                        }

                                        if(TDim>1)
                                        {
                                            double y = rNodePtr->GetCoordinate(1);
                                            if ((y >= 0.0   - Tol   && y <= 0.0   + Tol) ||
                                                (y >= rL[1] - Tol   && y <= rL[1] + Tol))
                                            {
                                                return true;
                                            }
                                        }
                                        if(TDim>2)
                                        {
                                            double z = rNodePtr->GetCoordinate(2);
                                            if ((z >=  0.0   - Tol   && z <= 0.0   + Tol) ||
                                                (z >=  rL[2] - Tol   && z <= rL[2] + Tol))
                                            {
                                                return true;
                                            }
                                        }
                                    }
                                    return false;
                                };  // GetBoundaryNodesLambda




    auto LambdaTimeDepBoundaryRH= [&MTCtrl,&tCtrl]
                                  (double rTime)->double
                                {
                                    if(rTime == 0.0)
                                    {
                                        return MTCtrl.InitialRelativeHumidity;
                                    }
                                    else
                                    {
                                        if(rTime< tCtrl.BC_TransitionTime)
                                        {
                                            return MTCtrl.InitialRelativeHumidity -
                                                   sin(rTime / tCtrl.BC_TransitionTime * 3.14 /2.0) *
                                                   (MTCtrl.InitialRelativeHumidity-MTCtrl.BoundaryEnvironmentalRH);
                                        }
                                        {
                                            return MTCtrl.BoundaryEnvironmentalRH;
                                        }
                                    }
                                };  //TimeDepBoundaryRHLambda


    SetupConstrainedNodeBoundaryElements<TDim>(S,
                                               LambdaGetBoundaryNodes,
                                               TI,
                                               LambdaTimeDepBoundaryRH);



    auto lambdaGetNodeLeftBottomFront = [rL](NuTo::NodeBase* rNodePtr) -> bool
                                {
                                    if(rNodePtr->GetNumDisplacements()==0)
                                        return false;
                                    double Tol = 1.e-6;
                                    if (rNodePtr->GetNumCoordinates()>0)
                                    {
                                        double x=0.0,
                                               y=0.0,
                                               z=0.0;
                                        x = rNodePtr->GetCoordinate(0);
                                        if(TDim>1)
                                            y = rNodePtr->GetCoordinate(1);
                                        if(TDim>2)
                                            z = rNodePtr->GetCoordinate(2);

                                        if (x >= 0.0   - Tol   && x <= 0.0   + Tol &&
                                            y >= 0.0   - Tol   && y <= 0.0   + Tol &&
                                            z >= 0.0   - Tol   && z <= 0.0   + Tol )
                                        {
                                            return true;
                                        }
                                    }
                                    return false;
                                };  // lambdaGetNodeLeftBottom

    auto lambdaGetNodeLeftTopFront = [rL](NuTo::NodeBase* rNodePtr) -> bool
                                {
                                    if(rNodePtr->GetNumDisplacements()==0)
                                        return false;
                                    double Tol = 1.e-6;
                                    if (rNodePtr->GetNumCoordinates()>0)
                                    {
                                        double x=0.0,
                                               y=0.0,
                                               z=0.0;
                                        x = rNodePtr->GetCoordinate(0);
                                        if(TDim>1)
                                            y = rNodePtr->GetCoordinate(1);
                                        if(TDim>2)
                                            z = rNodePtr->GetCoordinate(2);

                                        if (x >= 0.0   - Tol   && x <= 0.0   + Tol &&
                                            y >= 0.0   - Tol   && y <= 0.0   + Tol &&
                                            z >= rL[2] - Tol   && z <= rL[2] + Tol )
                                        {
                                            return true;
                                        }
                                    }
                                    return false;
                                };  // lambdaGetNodeLeftBottom

    auto lambdaGetNodeLeftBottomBack = [rL](NuTo::NodeBase* rNodePtr) -> bool
                                {
                                    if(rNodePtr->GetNumDisplacements()==0)
                                        return false;
                                    double Tol = 1.e-6;
                                    if (rNodePtr->GetNumCoordinates()>0)
                                    {
                                        double x=0.0,
                                               y=0.0,
                                               z=0.0;
                                        x = rNodePtr->GetCoordinate(0);
                                        if(TDim>1)
                                            y = rNodePtr->GetCoordinate(1);
                                        if(TDim>2)
                                            z = rNodePtr->GetCoordinate(2);

                                        if (x >= 0.0   - Tol   && x <= 0.0   + Tol &&
                                            y >= rL[1] - Tol   && y <= rL[1] + Tol &&
                                            z >= 0.0   - Tol   && z <= 0.0   + Tol )
                                        {
                                            return true;
                                        }
                                    }
                                    return false;
                                };  // lambdaGetNodeLeftBottom




    MeCtrl.AddConstraint<TDim>(TI,
                               lambdaGetNodeLeftBottomFront,
                               0);
    if(TDim>1)
    {
        MeCtrl.AddConstraint<TDim>(TI,
                                   lambdaGetNodeLeftBottomFront,
                                   1);

        MeCtrl.AddConstraint<TDim>(TI,
                                   lambdaGetNodeLeftBottomBack,
                                   0);
    }
    if(TDim>2)
    {
        MeCtrl.AddConstraint<TDim>(TI,
                                   lambdaGetNodeLeftBottomFront,
                                   2);

        MeCtrl.AddConstraint<TDim>(TI,
                                   lambdaGetNodeLeftBottomBack,
                                   2);

        MeCtrl.AddConstraint<TDim>(TI,
                                   lambdaGetNodeLeftTopFront,
                                   0);

        MeCtrl.AddConstraint<TDim>(TI,
                                   lambdaGetNodeLeftTopFront,
                                   1);
    }

    MTCtrl.SetupStaticData();
    S.NodeBuildGlobalDofs();

    SetupMultiProcessor(S);

    SetupVisualize(S);

    SetupTimeIntegration(TI,
                         tCtrl,
                         resultDir,
                         rStaggered);
    NuTo::Timer timer("shrinkagetest");
    TI.Solve(tCtrl.t_final);

    CheckMoistureTransportResults<TDim>(S,
                                        rN,
                                        rL);
    CheckMechanicsResults<TDim>(S);

}



//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//  Main
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

int main()
{
    std::map<NuTo::Node::eDof,NuTo::Interpolation::eTypeOrder> dofIPTMap;
    dofIPTMap[NuTo::Node::eDof::COORDINATES]            = NuTo::Interpolation::EQUIDISTANT1;
    dofIPTMap[NuTo::Node::eDof::DISPLACEMENTS]          = NuTo::Interpolation::EQUIDISTANT1;
    dofIPTMap[NuTo::Node::eDof::RELATIVEHUMIDITY]       = NuTo::Interpolation::EQUIDISTANT1;
    dofIPTMap[NuTo::Node::eDof::WATERVOLUMEFRACTION]    = NuTo::Interpolation::EQUIDISTANT1;


    ShrinkageTest<1>({3},
                     {0.08},
                     dofIPTMap);

    ShrinkageTest<2>({3,3},
                     {0.08,0.08},
                     dofIPTMap);

    ShrinkageTest<3>({3,3,3},
                     {0.08,0.08,0.08},
                     dofIPTMap);


    ShrinkageTest<1>({3},
                     {0.08},
                     dofIPTMap,
                     true);

    ShrinkageTest<2>({3,3},
                     {0.08,0.08},
                     dofIPTMap,
                     true);

    ShrinkageTest<3>({3,3,3},
                     {0.08,0.08,0.08},
                     dofIPTMap,
                     true);
}
