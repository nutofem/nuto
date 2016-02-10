
#include <nuto/base/Exception.h>
#include <nuto/mechanics/structures/unstructured/Structure.h>
//#include <nuto/mechanics/timeIntegration/CrankNicolson.h>
#include <nuto/mechanics/timeIntegration/CrankNicolsonEvaluate.h>

#include "boost/filesystem.hpp"

#include <iostream>
#include <sys/stat.h>
#include <sys/time.h>

#include <eigen3/Eigen/Core>


void TestConstraints()
{

}


int main()
{
    try
    {
        // %%%%%%%%%%%%%%%%%
        // Declare Variables
        // %%%%%%%%%%%%%%%%%

        double Height           = 0.2;
        double Length           = 2.0;

        double ElementHeight    = 0.025;
        double ElementLength    = 0.025;

        double Thickness        = 1.0;

        double Density          = 1.0;      // N/mm³
        double PoissonRatio     = 0.0;
        double YoungsModulus    = 30000.0;

        double SimulationTime   = 40.0;
        unsigned int TimeSteps  = 3;



        std::string     VTKFolder                       = "ResultStructureEvaluate";


        timeval         time_begin, time_end;


        // %%%%%%%%%%%%%%%%
        // Create structure
        // %%%%%%%%%%%%%%%%

        NuTo::Structure myStructure(2);

        myStructure.SetNumTimeDerivatives(2);
        myStructure.SetShowTime(false);





        // %%%%%%%%%%%%
        // Create nodes
        // %%%%%%%%%%%%

        unsigned int NumNodesX = static_cast<unsigned int>(floor(Length/ElementLength))+1;
        unsigned int NumNodesY = static_cast<unsigned int>(floor(Height/ElementHeight))+1;

        double DeltaX = Length/(NumNodesX-1);
        double DeltaY = Height/(NumNodesY-1);

        int NodeNum(0);
        NuTo::FullVector<double,Eigen::Dynamic> Coordinates(2);

        for (unsigned int iY=0; iY < NumNodesY; iY++)
        {
            for (unsigned int iX=0; iX < NumNodesX; iX++)
            {
                Coordinates(0) = iX * DeltaX;
                Coordinates(1) = iY * DeltaY;
                myStructure.NodeCreate(NodeNum,Coordinates);
                NodeNum++;
            }
        }





        // %%%%%%%%%%%%%%%%%%
        // Interpolation Type
        // %%%%%%%%%%%%%%%%%%

        int myInterpolationType = myStructure.InterpolationTypeCreate("Quad2D");
        myStructure.InterpolationTypeAdd(myInterpolationType, NuTo::Node::COORDINATES, NuTo::Interpolation::eTypeOrder::EQUIDISTANT1);
        myStructure.InterpolationTypeAdd(myInterpolationType, NuTo::Node::DISPLACEMENTS, NuTo::Interpolation::eTypeOrder::EQUIDISTANT1);





        // %%%%%%%%%%%%%%%
        // Create elements
        // %%%%%%%%%%%%%%%

        unsigned int NumElementsX = NumNodesX-1;
        unsigned int NumElementsY = NumNodesY-1;

        NuTo::FullVector<int,Eigen::Dynamic> Nodes(4);

        for (unsigned int iY=0; iY< NumElementsY; iY++)
        {
            for (unsigned int iX=0; iX < NumElementsX; iX++)
            {
                Nodes(0) = iX    +  iY    *NumNodesX;
                Nodes(1) = iX +1 +  iY    *NumNodesX;
                Nodes(2) = iX +1 + (iY+1) *NumNodesX;
                Nodes(3) = iX    + (iY+1) *NumNodesX;
                myStructure.ElementCreate(myInterpolationType, Nodes, std::string("ConstitutiveLawIp"), std::string("StaticData"));
            }
        }

        myStructure.ElementTotalConvertToInterpolationType();
        //myStructure.NodeBuildGlobalDofs();false

        myStructure.UseMaximumIndependentSets(true);
        myStructure.CalculateMaximumIndependentSets();





        // %%%%%%%%%%%%%%
        // Create Section
        // %%%%%%%%%%%%%%

        int mySection = myStructure.SectionCreate("Plane_Stress");
        myStructure.SectionSetThickness(mySection,Thickness);
        myStructure.ElementTotalSetSection(mySection);





        // %%%%%%%%%%%%%%%%%%%%%%%
        // Create Constitutive Law
        // %%%%%%%%%%%%%%%%%%%%%%%

        int myConstitutiveLaw = myStructure.ConstitutiveLawCreate("LinearElasticEngineeringStress");
        myStructure.ConstitutiveLawSetParameterDouble(myConstitutiveLaw,NuTo::Constitutive::eConstitutiveParameter::YOUNGS_MODULUS,YoungsModulus);
        myStructure.ConstitutiveLawSetParameterDouble(myConstitutiveLaw,NuTo::Constitutive::eConstitutiveParameter::POISSONS_RATIO,PoissonRatio);
        myStructure.ConstitutiveLawSetParameterDouble(myConstitutiveLaw,NuTo::Constitutive::eConstitutiveParameter::DENSITY,Density);
        myStructure.ElementTotalSetConstitutiveLaw(myConstitutiveLaw);





        // %%%%%%%%%%
        // Set groups
        // %%%%%%%%%%





        //left boundary
        int GRPNodes_Left = myStructure.GroupCreate("Nodes");
        int Direction = 0;
        double Min = 0. - 0.01 * ElementLength;
        double Max = 0. + 0.01 * ElementLength;;
        myStructure.GroupAddNodeCoordinateRange(GRPNodes_Left,Direction,Min,Max);

        //right boundary
        int GRPNodes_Right = myStructure.GroupCreate("Nodes");
        Direction = 0;
        Min = Length - 0.01 * ElementLength;
        Max = Length + 0.01 * ElementLength;
        myStructure.GroupAddNodeCoordinateRange(GRPNodes_Right,Direction,Min,Max);

        //allNodes
        int GrpNodes_All = myStructure.GroupCreate("Nodes");
        Direction=0;
        Min = 0.-0.01 * ElementLength;
        Max = Length + 0.01 * ElementLength;
        myStructure.GroupAddNodeCoordinateRange(GrpNodes_All,Direction,Min,Max);

        int GRPNodes_Center = myStructure.GroupCreate("Nodes");

        auto PickNodeFunction = [Length,Height](NuTo::NodeBase* rNodePtr) -> bool
                                {
                                    if (rNodePtr->GetNumCoordinates()>0)
                                    {
                                        double x = rNodePtr->GetCoordinate(0);
                                        double y = rNodePtr->GetCoordinate(1);
                                        if (x >= Length/2 -0.1 &&
                                            x <= Length/2 +0.1 &&
                                            y >= Height / 2.0)
                                        {
                                            return true;
                                        }
                                    }
                                    return false;
                                };

        myStructure.GroupAddNodeFunction(GRPNodes_Center,PickNodeFunction);



        int GRPNodes_BottomRight = myStructure.GroupCreate("Nodes");

        auto PickNodeFunction2 = [Length,Height](NuTo::NodeBase* rNodePtr) -> bool
                                {
                                    if (rNodePtr->GetNumCoordinates()>0)
                                    {
                                        double x = rNodePtr->GetCoordinate(0);
                                        double y = rNodePtr->GetCoordinate(1);
                                        if (x >= Length - Length/3 -0.1 &&
                                            x <= Length - Length/3 +0.1 &&
                                            y <= Height / 2.0)
                                        {
                                            return true;
                                        }
                                    }
                                    return false;
                                };

        myStructure.GroupAddNodeFunction(GRPNodes_BottomRight,PickNodeFunction2);










        // %%%%%%%%%%%%%%%
        // Set constraints
        // %%%%%%%%%%%%%%%


        //set a sinusoidal quarter wave
        double Period = 40.0;
        NuTo::FullMatrix<double,Eigen::Dynamic,Eigen::Dynamic> DispRHS(51,2);

        double t = 0.0;

        for (int i=0; i<DispRHS.GetNumRows()-1; i++)
        {


            t =  static_cast<double>(i) /
                 (static_cast<double>(DispRHS.GetNumRows())-2.) *
                 0.25 * Period;

            DispRHS(i,0) = t;
            DispRHS(i,1) = 0.1 * sin(t/Period * 2. * M_PI);

            //loadRHSFactor(count,0) = 0;
        }
        DispRHS(DispRHS.GetNumRows()-1,0) = DispRHS(DispRHS.GetNumRows()-2,0) + 1.0;
        DispRHS(DispRHS.GetNumRows()-1,1) = DispRHS(DispRHS.GetNumRows()-2,1);


        auto DispFunctionRHS =    [SimulationTime](double rTime)->double
                                  {
                                        double QuarterTime = SimulationTime/4.0;
                                        if(rTime>=2*QuarterTime)
                                        {
                                            if (rTime<=3*QuarterTime)
                                            {
                                                return -0.2 * (rTime-2*QuarterTime) / QuarterTime;
                                            }
                                            return -0.2;
                                        }
                                        return 0.0;

                                  };

        auto DispFunctionCenter = [SimulationTime](double rTime)->double
                                  {
                                        double QuarterTime = SimulationTime/4.0;
                                        double factor = 0.0;

                                        if (rTime<=1*QuarterTime)
                                        {
                                            factor = rTime/QuarterTime * 0.1;
                                        }
                                        else if(rTime<=3*QuarterTime)
                                        {
                                            factor = 0.1;
                                        }
                                        else
                                        {
                                            factor = 0.1 - (rTime-3*QuarterTime)/QuarterTime * 0.1;
                                        }

                                        return factor * sin(rTime);
                                  };




        NuTo::FullMatrix<double,Eigen::Dynamic,Eigen::Dynamic> DirectionX(2,1);
        DirectionX.SetValue(0,0,1.0);

        NuTo::FullMatrix<double,Eigen::Dynamic,Eigen::Dynamic> DirectionY(2,1);
        DirectionY.SetValue(1,0,1.0);

        NuTo::FullMatrix<double,Eigen::Dynamic,Eigen::Dynamic> DirectionXY(2,1);
        DirectionXY.SetValue(0,0,1.0);
        DirectionXY.SetValue(1,0,1.0);






        //myStructure.LoadCreateNodeGroupForce(0,GRPLeftBottomSupport,DirectionX, 1);

                                   myStructure.ConstraintLinearSetDisplacementNodeGroup(GRPNodes_Left           ,DirectionX,0);
                                   myStructure.ConstraintLinearSetDisplacementNodeGroup(GRPNodes_Left           ,DirectionY,0);
                                   myStructure.ConstraintLinearSetDisplacementNodeGroup(GRPNodes_Right          ,DirectionX,0);
        int ConstraintRightDisp  = myStructure.ConstraintLinearSetDisplacementNodeGroup(GRPNodes_Right          ,DirectionY,0);
        int ConstraintCenterDisp = myStructure.ConstraintLinearSetDisplacementNodeGroup(GRPNodes_Center         ,DirectionY,0);
        //                           myStructure.ConstraintLinearSetDisplacementNodeGroup(GRPNodes_BottomRight    ,DirectionY,0);





        // %%%%%%%%%%%%%
        // Visualization
        // %%%%%%%%%%%%%

#ifdef ENABLE_VISUALIZE
        //mkdir(VTKFolder.c_str(),0777);
        int visualizationGroup = myStructure.GroupCreate(NuTo::Groups::eGroupId::Elements);
        myStructure.GroupAddElementsTotal(visualizationGroup);

        myStructure.AddVisualizationComponent(visualizationGroup, NuTo::VisualizeBase::DISPLACEMENTS);
        myStructure.AddVisualizationComponent(visualizationGroup, NuTo::VisualizeBase::ENGINEERING_STRAIN);
        myStructure.AddVisualizationComponent(visualizationGroup, NuTo::VisualizeBase::ENGINEERING_STRESS);
        //myStructure.ExportVtkDataFileElements(VTKFolder+"/"+VTKFile,false);
#endif





        // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        // Create time integration scheme
        // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        NuTo::CrankNicolsonEvaluate myTimeIntegrationScheme(&myStructure);

        //myTimeIntegrationScheme.SetDampingCoefficientMass(0.05);
        //5 timesteps to capture the quarter wave (5/quarter wave)
        myTimeIntegrationScheme.SetMaxTimeStep(Period/20.);
        myTimeIntegrationScheme.SetTimeStep(SimulationTime/TimeSteps);
        //myTimeIntegrationScheme.SetMaxTimeStep(10);
        myTimeIntegrationScheme.SetMinTimeStep(0.001*myTimeIntegrationScheme.GetMaxTimeStep());


        myTimeIntegrationScheme.AddTimeDependentConstraint(ConstraintRightDisp,  DispFunctionRHS);
        myTimeIntegrationScheme.AddTimeDependentConstraint(ConstraintCenterDisp, DispFunctionCenter);
        //myTimeIntegrationScheme.AddTimeDependentConstraint(ConstraintRightDisp, DispRHS);
        //myTimeIntegrationScheme.SetTimeDependentConstraint(ConstraintRightDisp, DispRHS);



        //set output during the simulation to false
        myStructure.SetShowTime(false);
        myStructure.SetNumProcessors(4);

        myTimeIntegrationScheme.AddResultTime("Time");
        myTimeIntegrationScheme.AddResultGroupNodeForce("Forces_GroupNodes_Left" ,GRPNodes_Left);
        myTimeIntegrationScheme.AddResultGroupNodeForce("Forces_GroupNodes_Right",GRPNodes_Right);

        //set result directory
        bool deleteResultDirectoryFirst(true);
        myTimeIntegrationScheme.SetResultDirectory(VTKFolder,deleteResultDirectoryFirst);


        //solve (perform Newton raphson iteration
        gettimeofday(&time_begin, NULL);
        myTimeIntegrationScheme.Solve(SimulationTime);
        gettimeofday(&time_end, NULL);
        std::cout << "elapsed time : " << (time_end.tv_sec - time_begin.tv_sec) + (time_end.tv_usec - time_begin.tv_usec)/1000000.0<< " seconds" << std::endl;




        myStructure.Info();

    }
    catch (NuTo::Exception &e)
    {
        std::cout << e.what() << std::endl;
        return 1;
    }

    std::cout << "Test successful" << std::endl;
    return 0;
}
