#pragma once

#include "nuto/math/FullVector.h"
#include "nuto/mechanics/constitutive/ConstitutiveEnum.h"
#include "nuto/mechanics/elements/ElementBase.h"
#include "nuto/mechanics/elements/IpDataEnum.h"
#include "nuto/mechanics/nodes/NodeBase.h"
#include "nuto/mechanics/nodes/NodeEnum.h"
#include "nuto/mechanics/structures/unstructured/Structure.h"
#include "nuto/mechanics/integrationtypes/IntegrationTypeEnum.h"
#include "nuto/mechanics/timeIntegration/NewmarkDirect.h"
#include "nuto/mechanics/tools/MeshGenerator.h"
#include "nuto/mechanics/constitutive/staticData/DataMoistureTransport.h"
#include "nuto/mechanics/constitutive/laws/MoistureTransport.h"

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

#define RES_TOLERANCE 1e-18
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

/*---------------------------------------------*\
|*             moisture transport              *|
\*---------------------------------------------*/

struct MoistureTransportControl
{
    // general
    bool            EnableModiefiedTangentialStiffness  =   false;

    double          InitialRelativeHumidity                 =   1.0;    //VHIRTHAMTODO --- Set to 0.95
    double          InitialWaterVolumeFraction              =   0.;         //! --> Calculated from relative humidity
    double          MassExchangeRate                        =   3.42e-7;
    double          PoreVolumeFraction                           =   0.25;
    double          DiffusionCoefficientRH                  =   3.9e-12; //VHIRTHAMTODO ---> Vapor phase / Water phase ersetzen durch WVF und RH
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

    int             ConstitutiveLawID                       = 0;
    NuTo::Structure&    mS;

    // ctor
    MoistureTransportControl(NuTo::Structure& rS)
        :mS(rS)
    {
        // Values fitted from figure in Johannessons paper
        AdsorptionCoeffs(0) =  0.0;
        AdsorptionCoeffs(1) =  0.19692057340725558;
        AdsorptionCoeffs(2) = -0.28253538941816925;
        AdsorptionCoeffs(3) =  0.22661481601091368;

        DesorptionCoeffs(0) =  0.0;
        DesorptionCoeffs(1) =  0.26719233184420238;
        DesorptionCoeffs(2) = -0.41030868184510738;
        DesorptionCoeffs(3) =  0.32511635000090505;

        ConstitutiveLawID   = mS.ConstitutiveLawCreate("Moisture_Transport");
    }

public:
    void SetParametersConstitutiveLaw()
    {
        // set variables
        mS.ConstitutiveLawSetParameterBool    (ConstitutiveLawID  ,"ENABLE_MODIFIED_TANGENTIAL_STIFFNESS",EnableModiefiedTangentialStiffness);
        mS.ConstitutiveLawSetParameterBool    (ConstitutiveLawID  ,"enable_sorption_hysteresis",EnableSorptionHysteresis);

        mS.ConstitutiveLawSetParameterDouble  (ConstitutiveLawID  ,"BOUNDARY_DIFFUSION_COEFFICIENT_RH",BoundaryDiffusionCoefficientRH);
        mS.ConstitutiveLawSetParameterDouble  (ConstitutiveLawID  ,"BOUNDARY_DIFFUSION_COEFFICIENT_WV",BoundaryDiffusionCoefficientWV);
        mS.ConstitutiveLawSetParameterDouble  (ConstitutiveLawID  ,"DENSITY_WATER",DensityWater);
        mS.ConstitutiveLawSetParameterDouble  (ConstitutiveLawID  ,"DIFFUSION_COEFFICIENT_RH",DiffusionCoefficientRH);
        mS.ConstitutiveLawSetParameterDouble  (ConstitutiveLawID  ,"DIFFUSION_COEFFICIENT_WV",DiffusionCoefficientWV);
        mS.ConstitutiveLawSetParameterDouble  (ConstitutiveLawID  ,"DIFFUSION_EXPONENT_RH",DiffusionExponentRH);
        mS.ConstitutiveLawSetParameterDouble  (ConstitutiveLawID  ,"DIFFUSION_EXPONENT_WV",DiffusionExponentWV);
        mS.ConstitutiveLawSetParameterDouble  (ConstitutiveLawID  ,"GRADIENT_CORRECTION_ADSORPTION_DESORPTION",GradientCorrectionAdsorpionDesorption);
        mS.ConstitutiveLawSetParameterDouble  (ConstitutiveLawID  ,"GRADIENT_CORRECTION_DESORPTION_ADSORPTION",GradientCorrectionDesorptionAdsorpion);
        mS.ConstitutiveLawSetParameterDouble  (ConstitutiveLawID  ,"MASS_EXCHANGE_RATE",MassExchangeRate);
        mS.ConstitutiveLawSetParameterDouble  (ConstitutiveLawID  ,"PORE_VOLUME_FRACTION",PoreVolumeFraction);
        mS.ConstitutiveLawSetParameterDouble  (ConstitutiveLawID  ,"DENSITY_SATURATED_WATER_VAPOR",DensitySaturatedWaterVapor);

        mS.ConstitutiveLawSetParameterFullVectorDouble    (ConstitutiveLawID  ,"polynomial_COEFFICIENTS_ADSORPTION",AdsorptionCoeffs);
        mS.ConstitutiveLawSetParameterFullVectorDouble    (ConstitutiveLawID  ,"POLYNOMIAL_COEFFICIENTS_DESORPTION",DesorptionCoeffs);


        // Calculate equilibrium water volume fraction
        InitialWaterVolumeFraction   = mS.ConstitutiveLawGetEquilibriumWaterVolumeFraction(ConstitutiveLawID  ,
                                                                                           InitialRelativeHumidity,
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
                mS.NodeGetNodePtr(i)->Set(NuTo::Node::eDof::RELATIVEHUMIDITY, 0, InitialRelativeHumidity) ;
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
                auto& moistureData = mS.ElementGetElementPtr(i)->GetIPData()
                    .GetIPConstitutiveLaw(theIP)
                    .GetData<NuTo::MoistureTransport>()
                    .GetData();
                moistureData.SetLastSorptionCoeff(mS.ConstitutiveLawGetParameterFullVectorDouble(ConstitutiveLawID,NuTo::Constitutive::eConstitutiveParameter::POLYNOMIAL_COEFFICIENTS_DESORPTION));
                moistureData.SetCurrentSorptionCoeff(mS.ConstitutiveLawGetParameterFullVectorDouble(ConstitutiveLawID,NuTo::Constitutive::eConstitutiveParameter::POLYNOMIAL_COEFFICIENTS_DESORPTION));
                moistureData.SetLastRelHumValue(InitialRelativeHumidity);
                moistureData.SetDesorption(SorptionHistoryDesorption);
            }
        }
    }


    int SetupConstraint(int rBoundaryNodeID)
    {
        return  mS.ConstraintLinearSetRelativeHumidityNode(rBoundaryNodeID,1.0);
    }


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

    //rS.InterpolationTypeSetIntegrationType();

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
               std::array<double,TDim> rL);



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

//NuTo::MeshGenerator::MeshCylinder(rS,
//                                rSEC,
//                                rConsLaw,
//                                rIPT,
//                                {10,10,30},
//                                0.1,
//                                0.3);
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

template <int TDim>
inline void SetupStructure(NuTo::Structure& rS, std::string rTestName)
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
                                 const std::string& rResultDir)
{
    rTI.SetPerformLineSearch(false);
    rTI.SetVerboseLevel(0);
    rTI.SetToleranceResidual(NuTo::Node::eDof::RELATIVEHUMIDITY,RES_TOLERANCE);
    rTI.SetToleranceResidual(NuTo::Node::eDof::WATERVOLUMEFRACTION,RES_TOLERANCE);
    rTI.SetMaxNumIterations(MAX_ITERATION);

    rTI.SetTimeStep(rTC.delta_t);
    rTI.SetMinTimeStepPlot(rTC.t_write);

    rTI.SetResultDirectory(rResultDir,true);
}

/*---------------------------------------------*\
|*                 visualize                   *|
\*---------------------------------------------*/

inline void SetupVisualize(NuTo::Structure& rS)
{
#ifdef ENABLE_VISUALIZE
        int visGrp = rS.GroupCreate(NuTo::eGroupId::Elements);
        rS.GroupAddElementsTotal(visGrp);

        rS.AddVisualizationComponent(visGrp, NuTo::eVisualizeWhat::RELATIVE_HUMIDITY);
        rS.AddVisualizationComponent(visGrp, NuTo::eVisualizeWhat::WATER_VOLUME_FRACTION);
#endif // ENABLE_VISUALIZE
}
