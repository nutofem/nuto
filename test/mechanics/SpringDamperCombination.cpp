#include "nuto/math/FullMatrix.h"
#include "nuto/mechanics/constitutive/ConstitutiveEnum.h"
#include "nuto/mechanics/constitutive/laws/AdditiveOutput.h"
#include "nuto/mechanics/groups/GroupEnum.h"
#include "nuto/mechanics/interpolationtypes/InterpolationTypeEnum.h"
#include "nuto/mechanics/nodes/NodeBase.h"
#include "nuto/mechanics/nodes/NodeEnum.h"
#include "nuto/mechanics/structures/unstructured/Structure.h"
#include "nuto/mechanics/timeIntegration/TimeIntegrationBase.h"
#include "nuto/mechanics/timeIntegration/NewmarkBase.h"
#include "nuto/mechanics/timeIntegration/NewmarkDirect.h"
#include "nuto/mechanics/tools/MeshGenerator.h"
#include "nuto/visualize/VisualizeEnum.h"

#include <iostream>



//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//  Setup Preprocessor
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

/*---------------------------------------------*\
|*                  DEFINES                    *|
\*---------------------------------------------*/


// --- Test Parameters
// -------------------

#define SURFACELOAD 60.0e6

// --- Material Parameters
// -----------------------

#define LD_DAMPINGCOEFFICIENT 120.0e15
#define LE_YOUNGSMODULUS 30.0e9
#define LE_POISSONRATIO 0.2



// --- Processor/OpenMp
// --------------------

#ifdef _OPENMP
    #define TESTNUM_PROC 4
#elif HAVE_PARDISO
    #define TESTNUM_PROC 4
#else
    #define TESTNUM_PROC 1
#endif



// --- Time integration scheme
// ---------------------------
#define RES_TOLERANCE_MECHANICS 1e-6
#define MAX_ITERATION 20


//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//  Setup structs
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

/*---------------------------------------------*\
|*                    time                     *|
\*---------------------------------------------*/

struct TimeControl
{
    double          delta_t                         = 60.0;
    double          t_write                         = 60.0;
    double          t_final                         = 20.0 * 60.0;
};


//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//  Setup function
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



/*---------------------------------------------*\
|*                constraints                  *|
\*---------------------------------------------*/



template<int TDim>
int AddConstraint(NuTo::Structure& rS,
                  NuTo::NewmarkDirect& rTI,
                  std::function<bool(NuTo::NodeBase*)> rGetNodeFunction,
                  unsigned int rDirection,
                  double rValue = 0.0)
{
    assert(rValue <= TDim && "Direction isn't part of current dimension");
    int GRPNodesConstraint = rS.GroupCreate("Nodes");
    rS.GroupAddNodeFunction(GRPNodesConstraint,rGetNodeFunction);

    NuTo::FullMatrix<double,Eigen::Dynamic,Eigen::Dynamic> direction(TDim,1);
    direction.SetValue(rDirection, 0, 1.0);

    return rS.ConstraintLinearSetDisplacementNodeGroup(GRPNodesConstraint, direction,rValue);


}




/*---------------------------------------------*\
|*                   force                     *|
\*---------------------------------------------*/

template <int TDim>
void AddSurfaceLoad(NuTo::Structure& rS,
                    std::function<bool(NuTo::NodeBase*)> rGetNodeFunction)
{
    int GRPNodesSurfaceLoad = rS.GroupCreate("Nodes");
    rS.GroupAddNodeFunction(GRPNodesSurfaceLoad,rGetNodeFunction);

    NuTo::FullMatrix<double,Eigen::Dynamic,Eigen::Dynamic> direction(TDim,1);
    direction.SetValue(0, 0, 1.0);

    rS.LoadCreateNodeGroupForce(0,
                                GRPNodesSurfaceLoad,
                                direction,
                                SURFACELOAD);
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
                                 const std::string& rResultDir)
{
    rTI.SetPerformLineSearch(false);
    rTI.SetVerboseLevel(0);
    rTI.SetToleranceResidual(NuTo::Node::eDof::DISPLACEMENTS,RES_TOLERANCE_MECHANICS);
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
        rS.AddVisualizationComponent(visGrp, NuTo::eVisualizeWhat::DISPLACEMENTS);
        rS.AddVisualizationComponent(visGrp, NuTo::eVisualizeWhat::ENGINEERING_STRAIN);
        rS.AddVisualizationComponent(visGrp, NuTo::eVisualizeWhat::ENGINEERING_STRESS);
        rS.AddVisualizationComponent(visGrp, NuTo::eVisualizeWhat::PRINCIPAL_ENGINEERING_STRESS);
#endif // ENABLE_VISUALIZE
}

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//  Tests
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




template <int TDim>
void TestSpringDamperCombination(std::array<int,TDim> rN,
                                 std::array<double,TDim> rL,
                                 std::map<NuTo::Node::eDof,NuTo::Interpolation::eTypeOrder> rDofIPTMap)
{

    std::string testName = std::string("SpringDamperCombination") + std::to_string(TDim) +"D";

    std::cout << std::endl << "--------------------------------------------------------------------------"
              << std::endl << "Start test: "<< testName
              << std::endl << "--------------------------------------------------------------------------" << std::endl;




    // Allocate necessary resources
    NuTo::Structure S(TDim);
    NuTo::NewmarkDirect TI(&S);
    int CL_AO_ID = S.ConstitutiveLawCreate(NuTo::Constitutive::eConstitutiveType::ADDITIVE_OUTPUT);
    int CL_LD_ID = S.ConstitutiveLawCreate(NuTo::Constitutive::eConstitutiveType::LINEAR_DAMPING_ENGINEERING_STRESS);
    int CL_LE_ID = S.ConstitutiveLawCreate(NuTo::Constitutive::eConstitutiveType::LINEAR_ELASTIC_ENGINEERING_STRESS);

    NuTo::ConstitutiveBase* CL_AO_Ptr   = S.ConstitutiveLawGetConstitutiveLawPtr(CL_AO_ID);
    NuTo::ConstitutiveBase* CL_LD_Ptr   = S.ConstitutiveLawGetConstitutiveLawPtr(CL_LD_ID);
    NuTo::ConstitutiveBase* CL_LE_Ptr   = S.ConstitutiveLawGetConstitutiveLawPtr(CL_LE_ID);

    dynamic_cast<NuTo::AdditiveOutput*>(CL_AO_Ptr)->AddConstitutiveLaw(*CL_LE_Ptr);
    dynamic_cast<NuTo::AdditiveOutput*>(CL_AO_Ptr)->AddConstitutiveLaw(*CL_LD_Ptr);

    TimeControl tCtrl;
    tCtrl.t_final = 365.0 * 24.0 * 60.0 * 60.0;
    tCtrl.delta_t = tCtrl.t_final / 10.;
    tCtrl.t_write = tCtrl.delta_t;


    CL_LD_Ptr->SetParameterDouble(NuTo::Constitutive::eConstitutiveParameter::DAMPING_COEFFICIENT, LD_DAMPINGCOEFFICIENT);
    CL_LE_Ptr->SetParameterDouble(NuTo::Constitutive::eConstitutiveParameter::POISSONS_RATIO, LE_POISSONRATIO);
    CL_LE_Ptr->SetParameterDouble(NuTo::Constitutive::eConstitutiveParameter::YOUNGS_MODULUS, LE_YOUNGSMODULUS);

    SetupStructure(S,testName);
    int SEC = SetupSection<TDim>(S);
    int IPT = SetupInterpolationType<TDim>(S,rDofIPTMap);

    SetupMesh<TDim>(S,
                    SEC,
                    CL_AO_ID,
                    IPT,
                    rN,
                    rL);

    S.ElementTotalConvertToInterpolationType();
    S.NodeBuildGlobalDofs();


    // Add constraint on the leftern side
    auto lambdaGetNodesLeftSurface = [](NuTo::NodeBase* rNodePtr) -> bool
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
                                };  // lambdaGetNodeLeftBottom


    AddConstraint<TDim>(S,
                        TI,
                        lambdaGetNodesLeftSurface,
                        0);


    // Add force on the rightern side
    auto lambdaGetNodesRightSurface = [rL](NuTo::NodeBase* rNodePtr) -> bool
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
                                };  // lambdaGetNodeLeftBottom


    AddSurfaceLoad<TDim>(S,
                         lambdaGetNodesRightSurface);


    NuTo::FullMatrix<double,Eigen::Dynamic,Eigen::Dynamic> TimeDependentLoadFactor(3,2);
    TimeDependentLoadFactor(0,0) = 0.;
    TimeDependentLoadFactor(1,0) = 1.;
    TimeDependentLoadFactor(2,0) = 2.;
    TimeDependentLoadFactor(0,1) = 0.;
    TimeDependentLoadFactor(1,1) = 1.;
    TimeDependentLoadFactor(2,1) = 1.;

    TI.SetTimeDependentLoadCase(0,TimeDependentLoadFactor);

    SetupMultiProcessor(S);
    SetupVisualize(S);

    SetupTimeIntegration(TI,
                         tCtrl,
                         testName);

    TI.Solve(tCtrl.t_final);

    std::cout << "Test - PASSED!" << std::endl << std::endl;
}



//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//  Main
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

int main()
{
    std::map<NuTo::Node::eDof,NuTo::Interpolation::eTypeOrder> dofIPTMap = {
        {NuTo::Node::eDof::COORDINATES,   NuTo::Interpolation::eTypeOrder::EQUIDISTANT1},
        {NuTo::Node::eDof::DISPLACEMENTS, NuTo::Interpolation::eTypeOrder::EQUIDISTANT1}};


    TestSpringDamperCombination<1>( {3},
                                    {0.01},
                                    dofIPTMap);

    std::cout << "Everything went well -  NO errors." << std::endl;
    return 0;
}
