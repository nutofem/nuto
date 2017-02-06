#include "mechanics/constitutive/ConstitutiveEnum.h"
#include "mechanics/constitutive/laws/AdditiveOutput.h"
#include "mechanics/constitutive/laws/AdditiveInputImplicit.h"
#include "mechanics/groups/GroupEnum.h"
#include "mechanics/interpolationtypes/InterpolationTypeEnum.h"
#include "mechanics/nodes/NodeBase.h"
#include "mechanics/nodes/NodeEnum.h"
#include "mechanics/structures/unstructured/Structure.h"
#include "mechanics/timeIntegration/TimeIntegrationBase.h"
#include "mechanics/timeIntegration/NewmarkBase.h"
#include "mechanics/timeIntegration/NewmarkDirect.h"
#include "mechanics/mesh/MeshGenerator.h"
#include "visualize/VisualizeEnum.h"

#include <boost/foreach.hpp>

#include <iostream>
#include <math.h>


//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//  Setup Preprocessor
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

/*---------------------------------------------*\
|*                  DEFINES                    *|
\*---------------------------------------------*/


// --- Test Parameters
// -------------------

#define SURFACELOAD 60.0e6
#define DELTAT 30.0
#define TWRITE 60.0
#define TEND 40 * DELTAT


// --- Material Parameters
// -----------------------

#define LD_DAMPINGCOEFFICIENT 30.0e12
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
#define RES_TOLERANCE_MECHANICS 1e-5
#define MAX_ITERATION 20



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
    double          delta_t                         = DELTAT;
    double          t_write                         = TWRITE;
    double          t_final                         = TEND;
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

    return rS.ConstraintLinearSetDisplacementNodeGroup(GRPNodesConstraint, Eigen::Matrix<double, TDim, 1>::UnitX(), rValue);


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

    rS.LoadCreateNodeGroupForce(0,
                                GRPNodesSurfaceLoad,
                                Eigen::Matrix<double, TDim, 1>::UnitX(),
                                SURFACELOAD);
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
|*                nodal values                 *|
\*---------------------------------------------*/


template <int TDim>
void ApplyInitialNodalValues(NuTo::Structure& rS,
                             std::array<int,TDim> rN,
                             double rInitialStrainRate)
{
    unsigned int NNodes = rS.GetNumNodes();
    Eigen::VectorXd NodalStartValues(TDim);
    NodalStartValues.setZero();
    NodalStartValues(0) = rInitialStrainRate;

    for (unsigned int i=0; i<NNodes; i++)
    {
        if(rS.NodeGetNodePtr(i)->GetNum(NuTo::Node::eDof::DISPLACEMENTS) != 0)
        {
            NuTo::NodeBase* NodePtr =  rS.NodeGetNodePtr(i);
            Eigen::VectorXd NodeCoordinates = NodePtr->Get(NuTo::Node::eDof::COORDINATES);
            if(NodeCoordinates(0) == 0.0)
                continue;
            NodePtr->Set(NuTo::Node::eDof::DISPLACEMENTS,1,NodalStartValues * NodeCoordinates(0) );
        }
    }
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
                                 const int IDNodeRight)
{
    rTI.SetPerformLineSearch(false);
    rTI.SetVerboseLevel(0);
    rTI.SetToleranceResidual(NuTo::Node::eDof::DISPLACEMENTS,RES_TOLERANCE_MECHANICS);
    rTI.SetMaxNumIterations(MAX_ITERATION);

    rTI.SetTimeStep(rTC.delta_t);
    rTI.SetMinTimeStepPlot(rTC.t_write);

    rTI.SetResultDirectory(rResultDir,true);

    rTI.AddResultNodeDisplacements("Displacements",IDNodeRight);
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
//        rS.AddVisualizationComponent(visGrp, NuTo::eVisualizeWhat::ENGINEERING_STRAIN);
        rS.AddVisualizationComponent(visGrp, NuTo::eVisualizeWhat::ENGINEERING_STRESS);
        rS.AddVisualizationComponent(visGrp, NuTo::eVisualizeWhat::PRINCIPAL_ENGINEERING_STRESS);
#endif // ENABLE_VISUALIZE
}



//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//  Check Results
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



template<int TDim>
void CheckResultsSpringDamper(NuTo::Structure& rS,
                              std::array<double,TDim> rL)
{

    const NodeMap& nodePtrMap = rS.NodeGetNodeMap();
    BOOST_FOREACH(NodeMap::const_iterator::value_type it, nodePtrMap)
    {
        const NuTo::NodeBase* nodePtr = it.second;
        if(nodePtr->GetNum(NuTo::Node::eDof::DISPLACEMENTS)<1)
        {
            continue;   // Nodes without Displacements cant be checked
        }
        for(int i=0; i<TDim; ++i)
        {
            double coord   = nodePtr->Get(NuTo::Node::eDof::COORDINATES)[i];
            if(coord <= 0.)
                continue;
            double strain_numerical  = nodePtr->Get(NuTo::Node::eDof::DISPLACEMENTS)[i] / coord;
            double strain_theoretical = SURFACELOAD / LE_YOUNGSMODULUS * (1 - std::exp(-LE_YOUNGSMODULUS/LD_DAMPINGCOEFFICIENT * (TEND ))); // -delta_t/2.0:Because the first timestep produces an offset
                                                                                                                                                        // between theoretical solution and numerical solution which is delta_t/2
            double ErrorPercentage = std::abs(1-strain_numerical/strain_theoretical);
            const double tolerance = 5e-5;
            if(ErrorPercentage>tolerance)
                throw NuTo::Exception(__PRETTY_FUNCTION__,"Difference to theoretical solution is bigger than 0.005%!");
        }
    }
        std::cout << "Displacements correct!" << std::endl;
}


template<int TDim>
void CheckResultsSpringDamperSerial(NuTo::Structure& rS,
                              std::array<double,TDim> rL)
{

    const NodeMap& nodePtrMap = rS.NodeGetNodeMap();
    BOOST_FOREACH(NodeMap::const_iterator::value_type it, nodePtrMap)
    {
        const NuTo::NodeBase* nodePtr = it.second;
        if(nodePtr->GetNum(NuTo::Node::eDof::DISPLACEMENTS)<1)
        {
            continue;   // Nodes without Displacements cant be checked
        }
        for(int i=0; i<TDim; ++i)
        {
            double coord   = nodePtr->Get(NuTo::Node::eDof::COORDINATES)[i];
            if(coord <= 0.)
                continue;
            double strain_numerical  = nodePtr->Get(NuTo::Node::eDof::DISPLACEMENTS)[i] / coord;
            double strain_theoretical = SURFACELOAD / LE_YOUNGSMODULUS * (1 - std::exp(-LE_YOUNGSMODULUS/LD_DAMPINGCOEFFICIENT * (TEND - DELTAT/2.0))) +
                                        SURFACELOAD / (2* LE_YOUNGSMODULUS) * (1 - std::exp(-LE_YOUNGSMODULUS/LD_DAMPINGCOEFFICIENT * (TEND - DELTAT/2.0)));

            double ErrorPercentage = std::abs(1-strain_numerical/strain_theoretical);
            const double tolerance = 5e-3;
            if(ErrorPercentage>tolerance)
                throw NuTo::Exception(__PRETTY_FUNCTION__,"Difference to theoretical solution is bigger than 0.5%!");
        }
    }
        std::cout << "Displacements correct!" << std::endl;
}


//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//  Tests
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




template <int TDim>
void TestSpringDamperCombination(std::array<int,TDim> rN,
                                 std::array<double,TDim> rL)
{
    if(TDim>1)
        throw NuTo::Exception(__PRETTY_FUNCTION__,"2D and 3D are currently not supported!");

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

    CL_LD_Ptr->SetParameterDouble(NuTo::Constitutive::eConstitutiveParameter::DAMPING_COEFFICIENT, LD_DAMPINGCOEFFICIENT);
    CL_LE_Ptr->SetParameterDouble(NuTo::Constitutive::eConstitutiveParameter::POISSONS_RATIO, LE_POISSONRATIO);
    CL_LE_Ptr->SetParameterDouble(NuTo::Constitutive::eConstitutiveParameter::YOUNGS_MODULUS, LE_YOUNGSMODULUS);

    SetupStructure(S,testName);
    int SEC = SetupSection<TDim>(S);

    auto meshInfo = NuTo::MeshGenerator::Grid<TDim>(S, rL, rN);

    S.InterpolationTypeAdd(meshInfo.second, NuTo::Node::eDof::DISPLACEMENTS, NuTo::Interpolation::eTypeOrder::EQUIDISTANT1);

    S.ElementGroupSetSection(meshInfo.first, SEC);
    S.ElementGroupSetConstitutiveLaw(meshInfo.first, CL_AO_ID);

    S.ElementTotalConvertToInterpolationType();
    S.NodeBuildGlobalDofs();

//    ApplyInitialNodalValues<TDim>(S,
//                                  rN,
//                                  SURFACELOAD / LD_DAMPINGCOEFFICIENT);


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




//    NuTo::FullMatrix<double,Eigen::Dynamic,Eigen::Dynamic> TimeDependentLoadFactor(3,2);
//    TimeDependentLoadFactor(0,0) = 0.;
//    TimeDependentLoadFactor(1,0) = 1.;
//    TimeDependentLoadFactor(2,0) = 2.;
//    TimeDependentLoadFactor(0,1) = 1.;
//    TimeDependentLoadFactor(1,1) = 1.;
//    TimeDependentLoadFactor(2,1) = 1.;

//    TI.SetTimeDependentLoadCase(0,TimeDependentLoadFactor);

    SetupMultiProcessor(S);
    SetupVisualize(S);



    int IDNodeRight = -1;
    const auto& nodePtrMap = S.NodeGetNodeMap();
    BOOST_FOREACH(NodeMap::const_iterator::value_type it, nodePtrMap)
    {
        const NuTo::NodeBase* nodePtr = it.second;
        double coord   = nodePtr->Get(NuTo::Node::eDof::COORDINATES)[0];
        double tolerance = 1e-6;
        if(std::abs(coord-rL[0]) < tolerance)
        {
            IDNodeRight = S.NodeGetId(nodePtr);
            break;
        }
    }

    SetupTimeIntegration(TI,
                         tCtrl,
                         testName,
                         IDNodeRight);

    S.CalculateInitialValueRates(TI);
    TI.Solve(tCtrl.t_final);

    CheckResultsSpringDamper<TDim>(S,
                       rL);

    std::cout << "Test - PASSED!" << std::endl << std::endl;
}



template <int TDim>
void TestSpringDamperSerialChain(std::array<int,TDim> rN,
                                 std::array<double,TDim> rL)
{
    if(TDim>1)
        throw NuTo::Exception(__PRETTY_FUNCTION__,"2D and 3D are currently not supported!");

    std::string testName = std::string("SpringDamperSerialChain") + std::to_string(TDim) +"D";

    std::cout << std::endl << "--------------------------------------------------------------------------"
              << std::endl << "Start test: "<< testName
              << std::endl << "--------------------------------------------------------------------------" << std::endl;




    // Allocate necessary resources
    NuTo::Structure S(TDim);
    NuTo::NewmarkDirect TI(&S);

    int CL_AII_ID = S.ConstitutiveLawCreate(NuTo::Constitutive::eConstitutiveType::ADDITIVE_INPUT_IMPLICIT);
    int CL_AO1_ID = S.ConstitutiveLawCreate(NuTo::Constitutive::eConstitutiveType::ADDITIVE_OUTPUT);
    int CL_LD1_ID = S.ConstitutiveLawCreate(NuTo::Constitutive::eConstitutiveType::LINEAR_DAMPING_ENGINEERING_STRESS);
    int CL_LE1_ID = S.ConstitutiveLawCreate(NuTo::Constitutive::eConstitutiveType::LINEAR_ELASTIC_ENGINEERING_STRESS);
    int CL_AO2_ID = S.ConstitutiveLawCreate(NuTo::Constitutive::eConstitutiveType::ADDITIVE_OUTPUT);
    int CL_LD2_ID = S.ConstitutiveLawCreate(NuTo::Constitutive::eConstitutiveType::LINEAR_DAMPING_ENGINEERING_STRESS);
    int CL_LE2_ID = S.ConstitutiveLawCreate(NuTo::Constitutive::eConstitutiveType::LINEAR_ELASTIC_ENGINEERING_STRESS);


    NuTo::ConstitutiveBase* CL_AII_Ptr   = S.ConstitutiveLawGetConstitutiveLawPtr(CL_AII_ID);
    NuTo::ConstitutiveBase* CL_AO1_Ptr   = S.ConstitutiveLawGetConstitutiveLawPtr(CL_AO1_ID);
    NuTo::ConstitutiveBase* CL_LD1_Ptr   = S.ConstitutiveLawGetConstitutiveLawPtr(CL_LD1_ID);
    NuTo::ConstitutiveBase* CL_LE1_Ptr   = S.ConstitutiveLawGetConstitutiveLawPtr(CL_LE1_ID);
    NuTo::ConstitutiveBase* CL_AO2_Ptr   = S.ConstitutiveLawGetConstitutiveLawPtr(CL_AO2_ID);
    NuTo::ConstitutiveBase* CL_LD2_Ptr   = S.ConstitutiveLawGetConstitutiveLawPtr(CL_LD2_ID);
    NuTo::ConstitutiveBase* CL_LE2_Ptr   = S.ConstitutiveLawGetConstitutiveLawPtr(CL_LE2_ID);


    dynamic_cast<NuTo::AdditiveOutput*>(CL_AO1_Ptr)->AddConstitutiveLaw(*CL_LE1_Ptr);
    dynamic_cast<NuTo::AdditiveOutput*>(CL_AO1_Ptr)->AddConstitutiveLaw(*CL_LD1_Ptr);
    dynamic_cast<NuTo::AdditiveOutput*>(CL_AO2_Ptr)->AddConstitutiveLaw(*CL_LE2_Ptr);
    dynamic_cast<NuTo::AdditiveOutput*>(CL_AO2_Ptr)->AddConstitutiveLaw(*CL_LD2_Ptr);

    dynamic_cast<NuTo::AdditiveInputImplicit*>(CL_AII_Ptr)->AddConstitutiveLaw(*CL_AO1_Ptr);
    dynamic_cast<NuTo::AdditiveInputImplicit*>(CL_AII_Ptr)->AddConstitutiveLaw(*CL_AO2_Ptr);

    TimeControl tCtrl;

    CL_LD1_Ptr->SetParameterDouble(NuTo::Constitutive::eConstitutiveParameter::DAMPING_COEFFICIENT, LD_DAMPINGCOEFFICIENT);
    CL_LE1_Ptr->SetParameterDouble(NuTo::Constitutive::eConstitutiveParameter::POISSONS_RATIO, LE_POISSONRATIO);
    CL_LE1_Ptr->SetParameterDouble(NuTo::Constitutive::eConstitutiveParameter::YOUNGS_MODULUS, LE_YOUNGSMODULUS);

    CL_LD2_Ptr->SetParameterDouble(NuTo::Constitutive::eConstitutiveParameter::DAMPING_COEFFICIENT, LD_DAMPINGCOEFFICIENT * 2);
    CL_LE2_Ptr->SetParameterDouble(NuTo::Constitutive::eConstitutiveParameter::POISSONS_RATIO, LE_POISSONRATIO);
    CL_LE2_Ptr->SetParameterDouble(NuTo::Constitutive::eConstitutiveParameter::YOUNGS_MODULUS, LE_YOUNGSMODULUS * 2);

    SetupStructure(S,testName);
    int SEC = SetupSection<TDim>(S);

    auto meshInfo = NuTo::MeshGenerator::Grid<TDim>(S, rL, rN);

    S.InterpolationTypeAdd(meshInfo.second, NuTo::Node::eDof::DISPLACEMENTS, NuTo::Interpolation::eTypeOrder::EQUIDISTANT1);

    S.ElementGroupSetSection(meshInfo.first, SEC);
    S.ElementGroupSetConstitutiveLaw(meshInfo.first, CL_AII_ID);
    S.ElementTotalConvertToInterpolationType();
    S.NodeBuildGlobalDofs();
//    ApplyInitialNodalValues<TDim>(S,
//                                  rN,
//                                  SURFACELOAD / LD_DAMPINGCOEFFICIENT + SURFACELOAD / (2 * LD_DAMPINGCOEFFICIENT));
//    ApplyInitialNodalValues<TDim>(S,
//                                  rN,
//                                  0);


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




//    NuTo::FullMatrix<double,Eigen::Dynamic,Eigen::Dynamic> TimeDependentLoadFactor(3,2);
//    TimeDependentLoadFactor(0,0) = 0.;
//    TimeDependentLoadFactor(1,0) = 1.;
//    TimeDependentLoadFactor(2,0) = 2.;
//    TimeDependentLoadFactor(0,1) = 1.;
//    TimeDependentLoadFactor(1,1) = 1.;
//    TimeDependentLoadFactor(2,1) = 1.;

//    TI.SetTimeDependentLoadCase(0,TimeDependentLoadFactor);

    SetupMultiProcessor(S);
    SetupVisualize(S);



    int IDNodeRight = -1;
    const auto& nodePtrMap = S.NodeGetNodeMap();
    BOOST_FOREACH(NodeMap::const_iterator::value_type it, nodePtrMap)
    {
        const NuTo::NodeBase* nodePtr = it.second;
        double coord   = nodePtr->Get(NuTo::Node::eDof::COORDINATES)[0];
        double tolerance = 1e-6;
        if(std::abs(coord-rL[0]) < tolerance)
        {
            IDNodeRight = S.NodeGetId(nodePtr);
            break;
        }
    }

    SetupTimeIntegration(TI,
                         tCtrl,
                         testName,
                         IDNodeRight);


    S.CalculateInitialValueRates(TI);
    TI.Solve(tCtrl.t_final);

    CheckResultsSpringDamperSerial<TDim>(S,
                       rL);

    std::cout << "Test - PASSED!" << std::endl << std::endl;
}


//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//  Main
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

int main()
{
    TestSpringDamperCombination<1>( {10},
                                    {0.01});


    TestSpringDamperSerialChain<1>( {10},
                                    {0.01});

    std::cout << "Everything went well -  NO errors." << std::endl;
    return 0;
}
