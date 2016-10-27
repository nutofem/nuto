#include "nuto/mechanics/elements/IpDataEnum.h"
#include "nuto/mechanics/structures/unstructured/Structure.h"
#include "nuto/mechanics/timeIntegration/NewmarkDirect.h"
#include "nuto/mechanics/tools/MeshGenerator.h"
#include <array>
#include <boost/foreach.hpp>
#include "nuto/mechanics/constitutive/ConstitutiveEnum.h"
#include "nuto/mechanics/constitutive/laws/MoistureTransport.h"
#include "nuto/mechanics/constitutive/staticData/DataMoistureTransport.h"
#include "nuto/mechanics/constitutive/laws/AdditiveOutput.h"
#include "nuto/mechanics/constitutive/laws/AdditiveInputExplicit.h"
#include "nuto/mechanics/elements/ElementBase.h"
#include "nuto/mechanics/integrationtypes/IntegrationTypeEnum.h"
#include "nuto/mechanics/interpolationtypes/InterpolationTypeEnum.h"
#include "nuto/mechanics/nodes/NodeBase.h"
#include "nuto/mechanics/nodes/NodeEnum.h"
#ifdef ENABLE_VISUALIZE
#include "nuto/visualize/VisualizeEnum.h"
#include "nuto/mechanics/groups/GroupEnum.h"
#endif


//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//  Setup Preprocessor
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


/*---------------------------------------------*\
|*                  DEFINES                    *|
\*---------------------------------------------*/


// --- Time integration scheme
// ---------------------------
#define RES_TOLERANCE_MECHANICS 1e-4
#define RES_TOLERANCE_MOISTURE_TRANSPORT 1e-18
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

    double mYoungsModulus   =   30.e9;
    double mPoissonRatio    =   0.2;
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


            if(mS.NodeGetNodePtr(i)->GetNum(NuTo::Node::eDof::RELATIVEHUMIDITY) != 0)
            {
                mS.NodeGetNodePtr(i)->Set(NuTo::Node::eDof::RELATIVEHUMIDITY, 0,InitialRelativeHumidity) ;
            }
            if(mS.NodeGetNodePtr(i)->GetNum(NuTo::Node::eDof::WATERVOLUMEFRACTION) != 0)
            {
                mS.NodeGetNodePtr(i)->Set(NuTo::Node::eDof::WATERVOLUMEFRACTION, 0, InitialWaterVolumeFraction);
            }
        }
    }

    void SetupStaticData()
    {
        using namespace NuTo::Constitutive::StaticData;
        for (int i=0; i<mS.GetNumElements(); i++)
        {
            for (int theIP=0; theIP< mS.ElementGetElementPtr(i)->GetNumIntegrationPoints(); theIP++)
            {
                NuTo::ConstitutiveBase& lawbase
                    = mS.ElementGetElementPtr(i)->GetIPData().GetIPConstitutiveLaw(theIP).GetConstitutiveLaw();

                NuTo::AdditiveBase& lawAdditiveBase
                    = dynamic_cast<NuTo::AdditiveBase&>(lawbase);

                NuTo::Constitutive::IPConstitutiveLawBase& ipLawMoisture
                    = lawAdditiveBase.GetSublaw(1);

                NuTo::Constitutive::StaticData::DataMoistureTransport& moistureData
                    = ipLawMoisture.GetData<NuTo::MoistureTransport>().GetData();

                moistureData.SetLastSorptionCoeff(mMT.GetParameterFullVectorDouble(NuTo::Constitutive::eConstitutiveParameter::POLYNOMIAL_COEFFICIENTS_DESORPTION));
                moistureData.SetCurrentSorptionCoeff(mMT.GetParameterFullVectorDouble(NuTo::Constitutive::eConstitutiveParameter::POLYNOMIAL_COEFFICIENTS_DESORPTION));
                moistureData.SetLastRelHumValue(InitialRelativeHumidity);
                moistureData.SetDesorption(SorptionHistoryDesorption);
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
    controlNodeDofs.insert(NuTo::Node::eDof::RELATIVEHUMIDITY);

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
            elementPtr->SetIntegrationType(*rS.GetPtrIntegrationType(NuTo::eIntegrationType::IntegrationType0DBoundary));
            break;
        case 2:
            elementPtr->SetIntegrationType(*rS.GetPtrIntegrationType(NuTo::eIntegrationType::IntegrationType1D2NGauss2Ip));
            break;
        case 3:
            elementPtr->SetIntegrationType(*rS.GetPtrIntegrationType(NuTo::eIntegrationType::IntegrationType2D4NGauss4Ip));
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
        rS.InterpolationTypeSetIntegrationType(rIPT,NuTo::eIntegrationType::IntegrationType1D2NGauss2Ip);
        break;
    case 2:
        rS.InterpolationTypeSetIntegrationType(rIPT,NuTo::eIntegrationType::IntegrationType2D4NGauss4Ip);
        break;
    case 3:
        rS.InterpolationTypeSetIntegrationType(rIPT,NuTo::eIntegrationType::IntegrationType3D8NGauss2x2x2Ip);
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
        int visGrp = rS.GroupCreate(NuTo::eGroupId::Elements);
        rS.GroupAddElementsTotal(visGrp);
        rS.AddVisualizationComponent(visGrp, NuTo::eVisualizeWhat::DISPLACEMENTS);
        rS.AddVisualizationComponent(visGrp, NuTo::eVisualizeWhat::RELATIVE_HUMIDITY);
        rS.AddVisualizationComponent(visGrp, NuTo::eVisualizeWhat::WATER_VOLUME_FRACTION);
        rS.AddVisualizationComponent(visGrp, NuTo::eVisualizeWhat::PRINCIPAL_ENGINEERING_STRESS);
#endif // ENABLE_VISUALIZE
}








//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//  Check Results
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


void CheckMechanicsResults(NuTo::Structure& rS)
{
    const NodeMap& nodePtrMap = rS.NodeGetNodeMap();
    BOOST_FOREACH(NodeMap::const_iterator::value_type it, nodePtrMap)
    {
        const NuTo::NodeBase* nodePtr = it.second;
        if(nodePtr->GetNum(NuTo::Node::eDof::DISPLACEMENTS)<1)
        {
            continue;   // Nodes without Displacements cant be checked
        }
        double coordX = nodePtr->Get(NuTo::Node::eDof::COORDINATES)[0];
        double dispX  = nodePtr->Get(NuTo::Node::eDof::DISPLACEMENTS)[0];
        double diff   = std::abs(coordX * 0.1) - std::abs(dispX);
        constexpr const double tolerance = 1e-6;
        if((diff >tolerance || diff < -tolerance))// && coordX > 0)
        {
            throw NuTo::Exception(__PRETTY_FUNCTION__,"One ore more calculated Displacement is not correct");
        }
    }
        std::cout << "Displacements correct!" << std::endl;
}

template<int TDim>
void CheckMoistureTransportResults(NuTo::Structure& rS,
                                   std::array<int,TDim> rN,
                                   std::array<double,TDim> rL)
{
    assert((rN[0] == 16)   && "Only 16 elements in x-direction allowed for this test --- needed for direct comparison to values given in Johannesson 2010");
    assert((rL[0] == 0.16) && "The length in flow direction (x) must be 0.16m for direct comparison with paper values");

    // values fitted from Johannesson and Nyman(2010)
    NuTo::FullVector<double,Eigen::Dynamic> PaperValues(17);
    PaperValues(0)  = 0.06;
    PaperValues(1)  = 0.097;
    PaperValues(2)  = 0.116;
    PaperValues(3)  = 0.129;
    PaperValues(4)  = 0.138;
    PaperValues(5)  = 0.146;
    PaperValues(6)  = 0.148;
    PaperValues(7)  = 0.151;
    PaperValues(8)  = 0.152;
    PaperValues(9)  = PaperValues(7);
    PaperValues(10) = PaperValues(6);
    PaperValues(11) = PaperValues(5);
    PaperValues(12) = PaperValues(4);
    PaperValues(13) = PaperValues(3);
    PaperValues(14) = PaperValues(2);
    PaperValues(15) = PaperValues(1);
    PaperValues(16) = PaperValues(0);



    constexpr const double tolerance = 0.005; // Tolerance because not all necessary value (sorption curve) are given in the paper and must be approximated
    double deltaL = rL[0]/rN[0];
    unsigned int numMismatchingValues = 0;


    const NodeMap& nodePtrMap = rS.NodeGetNodeMap();
    BOOST_FOREACH(NodeMap::const_iterator::value_type it, nodePtrMap)
    {
        const NuTo::NodeBase* nodePtr = it.second;


        if(nodePtr->GetNum(NuTo::Node::eDof::WATERVOLUMEFRACTION)<1)
        {
            continue;   // Nodes without WVF cant be checked --- for example boundary control node
        }

        double relevantNodeCoord = nodePtr->Get(NuTo::Node::eDof::COORDINATES)[0];
        int relevantIndex = static_cast<int>(std::round(relevantNodeCoord / deltaL));

        double nodalWVF = nodePtr->Get(NuTo::Node::eDof::WATERVOLUMEFRACTION)[0];
        double paperWVF = PaperValues[relevantIndex];
        if(std::abs(nodalWVF-paperWVF) > tolerance)
        {
            std::cout << "nodalWVF nuto:  " << nodalWVF << std::endl;
            std::cout << "nodalWVF paper: " << paperWVF << std::endl;
            ++numMismatchingValues;
        }
    }
    if(numMismatchingValues>0)
    {
        throw NuTo::Exception(__PRETTY_FUNCTION__,"One or more calculated water volume fraction values exceed the tolerance when compared to reference values.");
    }
    std::cout << "Water volume fraction correct!" << std::endl;
}



//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//  Tests
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



//! @brief performs a simulation in the desired dimension
//! @param rN: array with number of elements in each direction
//! @param rL: array with length of elements in each direction
template<int TDim>
void AdditiveOutputTest(std::array<int,TDim> rN,
                        std::array<double,TDim> rL,
                        std::map<NuTo::Node::eDof,NuTo::Interpolation::eTypeOrder> rDofIPTMap,
                        bool rStaggered = false)
{
    assert((rN[0] == 16)   && "Only 16 elements in x-direction allowed for this test --- needed for direct comparison to values given in Johannesson 2010");
    assert((rL[0] == 0.16) && "The length in flow direction (x) must be 0.16m for direct comparison with paper values");


    std::string testName = std::string("AdditiveOutput") + std::to_string(TDim) +"D";
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
    int CL_MT_ID   = S.ConstitutiveLawCreate(NuTo::Constitutive::eConstitutiveType::MOISTURE_TRANSPORT);
    int CL_CLAL_ID = S.ConstitutiveLawCreate(NuTo::Constitutive::eConstitutiveType::ADDITIVE_OUTPUT);



    NuTo::ConstitutiveBase* CL_LE_Ptr   = S.ConstitutiveLawGetConstitutiveLawPtr(CL_LE_ID);
    NuTo::ConstitutiveBase* CL_MT_Ptr   = S.ConstitutiveLawGetConstitutiveLawPtr(CL_MT_ID);
    NuTo::AdditiveOutput* CL_CLAL_Ptr =
        static_cast<NuTo::AdditiveOutput*>(S.ConstitutiveLawGetConstitutiveLawPtr(CL_CLAL_ID));

    TimeControl tCtrl;
    tCtrl.t_final = 293.0 * 24.0 * 60.0 * 60.0;
    tCtrl.t_write = tCtrl.t_final;
    tCtrl.delta_t = tCtrl.t_final / 5.0;


    MoistureTransportControl MTCtrl(S,*CL_MT_Ptr);
    MTCtrl.InitialRelativeHumidity = 0.95;
    MTCtrl.MassExchangeRate = 3.42e-7;
    MTCtrl.DiffusionCoefficientRH = 3.9e-10;
    MTCtrl.BoundaryEnvironmentalRH = 0.40;
    MTCtrl.SetParametersConstitutiveLaw();

    MechanicsControl MeCtrl(S,*CL_LE_Ptr);
    MeCtrl.SetParametersConstitutiveLaw();

    CL_CLAL_Ptr->AddConstitutiveLaw(*CL_LE_Ptr);
    CL_CLAL_Ptr->AddConstitutiveLaw(*CL_MT_Ptr);


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
                                    if (rNodePtr->GetNum(NuTo::Node::eDof::COORDINATES)>0)
                                    {
                                        double x = rNodePtr->Get(NuTo::Node::eDof::COORDINATES)[0];
                                        if ((x >= 0.0   - Tol   && x <= 0.0   + Tol) ||
                                            (x >= rL[0] - Tol   && x <= rL[0] + Tol))
                                        {
                                            return true;
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




    auto lambdaGetNodesLeftSide = [rL](NuTo::NodeBase* rNodePtr) -> bool
                                {
                                    if(rNodePtr->GetNum(NuTo::Node::eDof::DISPLACEMENTS)==0)
                                        return false;
                                    double Tol = 1.e-6;
                                    if (rNodePtr->GetNum(NuTo::Node::eDof::COORDINATES)>0)
                                    {
                                        double x = rNodePtr->Get(NuTo::Node::eDof::COORDINATES)[0];
                                        if (x >= 0.0   - Tol   && x <= 0.0   + Tol)
                                        {
                                            return true;
                                        }
                                    }
                                    return false;
                                };  // lambdaGetNodesLeftSide

    auto lambdaGetNodeLeftBottom = [rL](NuTo::NodeBase* rNodePtr) -> bool
                                {
                                    if(rNodePtr->GetNum(NuTo::Node::eDof::DISPLACEMENTS)==0)
                                        return false;
                                    double Tol = 1.e-6;
                                    if (rNodePtr->GetNum(NuTo::Node::eDof::COORDINATES)>0)
                                    {
                                        double x=0.0,
                                               y=0.0,
                                               z=0.0;
                                        x = rNodePtr->Get(NuTo::Node::eDof::COORDINATES)[0];
                                        if(TDim>1)
                                            y = rNodePtr->Get(NuTo::Node::eDof::COORDINATES)[1];
                                        if(TDim>2)
                                            z = rNodePtr->Get(NuTo::Node::eDof::COORDINATES)[2];

                                        if (x >= 0.0   - Tol   && x <= 0.0   + Tol &&
                                            y >= 0.0   - Tol   && y <= 0.0   + Tol &&
                                            z >= 0.0   - Tol   && z <= 0.0   + Tol )
                                        {
                                            return true;
                                        }
                                    }
                                    return false;
                                };  // lambdaGetNodeLeftBottom

    auto lambdaGetNodesRightSide = [rL](NuTo::NodeBase* rNodePtr) -> bool
                                {
                                    if(rNodePtr->GetNum(NuTo::Node::eDof::DISPLACEMENTS)==0)
                                        return false;
                                    double Tol = 1.e-6;
                                    if (rNodePtr->GetNum(NuTo::Node::eDof::COORDINATES)>0)
                                    {
                                        double x = rNodePtr->Get(NuTo::Node::eDof::COORDINATES)[0];
                                        if (x >= rL[0] - Tol   && x <= rL[0] + Tol)
                                        {
                                            return true;
                                        }
                                    }
                                    return false;
                                };  // lambdaGetNodesRightSide


    auto lambdaDisplacementRightSide = [tCtrl,rL](double rTime)->double
                              {
                                    return -0.1 * rL[0] * rTime / tCtrl.t_final;
                              }; // lambdaDisplacementRightSide



    MeCtrl.AddConstraint<TDim>(TI,
                               lambdaGetNodesLeftSide,
                               0);
    if(TDim>1)
        MeCtrl.AddConstraint<TDim>(TI,
                                   lambdaGetNodeLeftBottom,
                                   1);
    if(TDim>2)
        MeCtrl.AddConstraint<TDim>(TI,
                                   lambdaGetNodeLeftBottom,
                                   2);
    MeCtrl.AddTimeDependentConstraint<TDim>(TI,
                                            lambdaGetNodesRightSide,
                                            lambdaDisplacementRightSide,
                                            0);
    MTCtrl.SetupStaticData();
    S.NodeBuildGlobalDofs();

    SetupMultiProcessor(S);

    SetupVisualize(S);

    SetupTimeIntegration(TI,
                         tCtrl,
                         resultDir,
                         rStaggered);
    TI.Solve(tCtrl.t_final);

    CheckMechanicsResults(S);
    CheckMoistureTransportResults<TDim>(S,
                                        rN,
                                        rL);

}








template<int TDim>
void AdditiveInputImplicitTest(std::array<int,TDim> rN,
                               std::array<double,TDim> rL,
                               std::map<NuTo::Node::eDof,NuTo::Interpolation::eTypeOrder> rDofIPTMap,
                               bool rStaggered = false)
{


    std::string testName = std::string("ConstitutiveLawsAdditiveInputImplicit") + std::to_string(TDim) +"D";
    if(rStaggered)
        testName += "_staggered";
    std::string resultDir = std::string("./MultipleConstitutiveLaws_") + testName;

    std::cout << std::endl << "--------------------------------------------------------------------------"
              << std::endl << "Start test: "<< testName
              << std::endl << "--------------------------------------------------------------------------" << std::endl;

    // Allocate neccessary stuff
    NuTo::Structure S(TDim);
    NuTo::NewmarkDirect TI(&S);
    int CL_LE1_ID   = S.ConstitutiveLawCreate(NuTo::Constitutive::eConstitutiveType::LINEAR_ELASTIC_ENGINEERING_STRESS);
    int CL_LE2_ID   = S.ConstitutiveLawCreate(NuTo::Constitutive::eConstitutiveType::LINEAR_ELASTIC_ENGINEERING_STRESS);
    int CL_AII_ID   = S.ConstitutiveLawCreate(NuTo::Constitutive::eConstitutiveType::ADDITIVE_INPUT_IMPLICIT);


    NuTo::ConstitutiveBase* CL_LE1_Ptr = S.ConstitutiveLawGetConstitutiveLawPtr(CL_LE1_ID);
    NuTo::ConstitutiveBase* CL_LE2_Ptr = S.ConstitutiveLawGetConstitutiveLawPtr(CL_LE2_ID);
    NuTo::AdditiveInputExplicit* CL_AII_Ptr =
        static_cast<NuTo::AdditiveInputExplicit*>(S.ConstitutiveLawGetConstitutiveLawPtr(CL_AII_ID));

    CL_AII_Ptr->AddConstitutiveLaw(*CL_LE1_Ptr);
    CL_AII_Ptr->AddConstitutiveLaw(*CL_LE2_Ptr);

    MechanicsControl MeCtrl1(S,*CL_LE1_Ptr);
    MechanicsControl MeCtrl2(S,*CL_LE2_Ptr);
    MeCtrl1.mYoungsModulus = 30e9;
    MeCtrl2.mYoungsModulus = 30e9;
    MeCtrl1.SetParametersConstitutiveLaw();
    MeCtrl2.SetParametersConstitutiveLaw();

    TimeControl tCtrl;
    tCtrl.t_final = 293.0 * 24.0 * 60.0 * 60.0;
    tCtrl.t_write = tCtrl.t_final;
    tCtrl.delta_t = tCtrl.t_final / 5.0;

    SetupStructure(S,testName);
    int SEC = SetupSection<TDim>(S);
    int IPT = SetupInterpolationType<TDim>(S,rDofIPTMap);

    SetupMesh<TDim>(S,
                    SEC,
                    CL_AII_ID,
                    IPT,
                    rN,
                    rL);


    S.ElementTotalConvertToInterpolationType(); //old used values 1.0e-12,0.001


    auto lambdaGetNodesLeftSide = [](NuTo::NodeBase* rNodePtr) -> bool
                                {
                                    if(rNodePtr->GetNum(NuTo::Node::eDof::DISPLACEMENTS)==0)
                                        return false;
                                    double Tol = 1.e-6;
                                    if (rNodePtr->GetNum(NuTo::Node::eDof::COORDINATES)>0)
                                    {
                                        double x = rNodePtr->Get(NuTo::Node::eDof::COORDINATES)[0];
                                        if (x >= 0.0   - Tol   && x <= 0.0   + Tol)
                                        {
                                            return true;
                                        }
                                    }
                                    return false;
                                };  // lambdaGetNodesLeftSide


    auto lambdaGetNodesRightSide = [&rL](NuTo::NodeBase* rNodePtr) -> bool
                                {
                                    if(rNodePtr->GetNum(NuTo::Node::eDof::DISPLACEMENTS)==0)
                                        return false;
                                    double Tol = 1.e-6;
                                    if (rNodePtr->GetNum(NuTo::Node::eDof::COORDINATES)>0)
                                    {
                                        double x = rNodePtr->Get(NuTo::Node::eDof::COORDINATES)[0];
                                        if (x >= rL[0]   - Tol   && x <= rL[0]   + Tol)
                                        {
                                            return true;
                                        }
                                    }
                                    return false;
                                };  // lambdaGetNodesLeftSide


    auto lambdaGetNodeLeftBottom = [](NuTo::NodeBase* rNodePtr) -> bool
                                {
                                    if(rNodePtr->GetNum(NuTo::Node::eDof::DISPLACEMENTS)==0)
                                        return false;
                                    double Tol = 1.e-6;
                                    if (rNodePtr->GetNum(NuTo::Node::eDof::COORDINATES)>0)
                                    {
                                        double x=0.0,
                                               y=0.0,
                                               z=0.0;
                                        x = rNodePtr->Get(NuTo::Node::eDof::COORDINATES)[0];
                                        if(TDim>1)
                                            y = rNodePtr->Get(NuTo::Node::eDof::COORDINATES)[1];
                                        if(TDim>2)
                                            z = rNodePtr->Get(NuTo::Node::eDof::COORDINATES)[2];

                                        if (x >= 0.0   - Tol   && x <= 0.0   + Tol &&
                                            y >= 0.0   - Tol   && y <= 0.0   + Tol &&
                                            z >= 0.0   - Tol   && z <= 0.0   + Tol )
                                        {
                                            return true;
                                        }
                                    }
                                    return false;
                                };  // lambdaGetNodeLeftBottom




    auto lambdaDisplacementRightSide = [tCtrl,rL](double rTime)->double
                              {
                                    return -0.1 * rL[0] * rTime / tCtrl.t_final;
                              }; // lambdaDisplacementRightSide



    MeCtrl1.AddConstraint<TDim>(TI,
                               lambdaGetNodesLeftSide,
                               0);
    if(TDim>1)
        MeCtrl1.AddConstraint<TDim>(TI,
                                   lambdaGetNodeLeftBottom,
                                   1);
    if(TDim>2)
        MeCtrl1.AddConstraint<TDim>(TI,
                                   lambdaGetNodeLeftBottom,
                                   2);

    MeCtrl1.AddTimeDependentConstraint<TDim>(TI,
                                            lambdaGetNodesRightSide,
                                            lambdaDisplacementRightSide,
                                            0);

    S.NodeBuildGlobalDofs(); //<--- possible memory leak!!! (Valgrind!!!)

    SetupMultiProcessor(S);

#ifdef ENABLE_VISUALIZE
        int visGrp = S.GroupCreate(NuTo::eGroupId::Elements);
        S.GroupAddElementsTotal(visGrp);
        S.AddVisualizationComponent(visGrp, NuTo::eVisualizeWhat::DISPLACEMENTS);
        S.AddVisualizationComponent(visGrp, NuTo::eVisualizeWhat::PRINCIPAL_ENGINEERING_STRESS);
#endif // ENABLE_VISUALIZE

    SetupTimeIntegration(TI,
                         tCtrl,
                         resultDir,
                         rStaggered);
    TI.Solve(tCtrl.t_final);

}



//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//  Main
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

int main()
{
    std::map<NuTo::Node::eDof,NuTo::Interpolation::eTypeOrder> dofIPTMap;
    dofIPTMap[NuTo::Node::eDof::COORDINATES]            = NuTo::Interpolation::eTypeOrder::EQUIDISTANT1;
    dofIPTMap[NuTo::Node::eDof::DISPLACEMENTS]          = NuTo::Interpolation::eTypeOrder::EQUIDISTANT1;
    dofIPTMap[NuTo::Node::eDof::RELATIVEHUMIDITY]       = NuTo::Interpolation::eTypeOrder::EQUIDISTANT1;
    dofIPTMap[NuTo::Node::eDof::WATERVOLUMEFRACTION]    = NuTo::Interpolation::eTypeOrder::EQUIDISTANT1;

    AdditiveOutputTest<1>({16},
                          {0.16},
                          dofIPTMap);

    AdditiveOutputTest<2>({16,2},
                          {0.16,0.02},
                          dofIPTMap);

    AdditiveOutputTest<3>({16,2,2},
                          {0.16,0.02,0.02},
                          dofIPTMap);

    AdditiveOutputTest<1>({16},
                          {0.16},
                          dofIPTMap,
                          true);

    AdditiveOutputTest<2>({16,2},
                          {0.16,0.02},
                          dofIPTMap,
                          true);

    AdditiveOutputTest<3>({16,2,2},
                          {0.16,0.02,0.02},
                          dofIPTMap,
                          true);


// Solver (MUMPS / PARDISO aren't thread save! Find other solution in constitutive law to solve local system)
#ifndef _OPENMP
    dofIPTMap.clear();
    dofIPTMap[NuTo::Node::eDof::COORDINATES]            = NuTo::Interpolation::eTypeOrder::EQUIDISTANT1;
    dofIPTMap[NuTo::Node::eDof::DISPLACEMENTS]          = NuTo::Interpolation::eTypeOrder::EQUIDISTANT1;

    AdditiveInputImplicitTest<1>({16},
                                 {0.16},
                                 dofIPTMap);

    AdditiveInputImplicitTest<2>({16,2},
                                 {0.16,0.02},
                                 dofIPTMap);

    AdditiveInputImplicitTest<3>({16,2,2},
                                 {0.16,0.02,0.02},
                                 dofIPTMap);
#endif
}
