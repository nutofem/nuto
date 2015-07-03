
#include <nuto/base/Exception.h>
#include <nuto/mechanics/structures/unstructured/Structure.h>
#include <nuto/mechanics/timeIntegration/NewmarkDirect.h>

#include "boost/filesystem.hpp"

#include <iostream>
#include <sys/stat.h>


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

        double ElementHeight    = 0.1;
        double ElementLength    = 0.1;

        double Thickness        = 1.0;

        double Density          = 1.0;      // N/mm³
        double PoissonRatio     = 0.0;
        double YoungsModulus    = 30000.0;

        double SimulationTime   = 10.0;
        unsigned int TimeSteps  = 50;



        std::string     VTKFile                         = "ResultsStructureEvaluate.vtk";
        std::string     VTKFolder                       = "Result";





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
        //myStructure.NodeBuildGlobalDofs();

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
        myStructure.ConstitutiveLawSetYoungsModulus(myConstitutiveLaw,YoungsModulus);
        myStructure.ConstitutiveLawSetPoissonsRatio(myConstitutiveLaw,PoissonRatio);
        myStructure.ConstitutiveLawSetDensity(myConstitutiveLaw,Density);
        myStructure.ElementTotalSetConstitutiveLaw(myConstitutiveLaw);





        // %%%%%%%%%%
        // Set groups
        // %%%%%%%%%%

        // bottom boundary
        int GRPNodes_Bottom = myStructure.GroupCreate("Nodes");
        int Direction=1;
        double Min= 0. - 0.01 * ElementHeight;
        double Max= 0. + 0.01 * ElementHeight;
        myStructure.GroupAddNodeCoordinateRange(GRPNodes_Bottom,Direction,Min,Max);

        // top boundary
        int GRPNodes_Top = myStructure.GroupCreate("Nodes");
        Direction=1;
        Min = Height - 0.01 * ElementHeight;
        Max = Height + 0.01 * ElementHeight;
        myStructure.GroupAddNodeCoordinateRange(GRPNodes_Top,Direction,Min,Max);

        //left boundary
        int GRPNodes_Left = myStructure.GroupCreate("Nodes");
        Direction = 0;
        Min = 0. - 0.01 * ElementLength;
        Max = 0. + 0.01 * ElementLength;;
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

        //intersect the groups for the left bottom node
        int GRPLeftBottomSupport = myStructure.GroupIntersection(GRPNodes_Bottom,GRPNodes_Left);
        if (myStructure.GroupGetNumMembers(GRPLeftBottomSupport)!=1)
        {
            std::cout << "group for bottom left boundary node has " << myStructure.GroupGetNumMembers(GRPLeftBottomSupport) << " members."<< std::endl;
            exit(-1);
        }

        int GRPLeftTopSupport = myStructure.GroupIntersection(GRPNodes_Top,GRPNodes_Left);
        if (myStructure.GroupGetNumMembers(GRPLeftTopSupport)!=1)
        {
            std::cout << "group for top left boundary node has " << myStructure.GroupGetNumMembers(GRPLeftTopSupport) << " members."<< std::endl;
            exit(-1);
        }







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



        NuTo::FullMatrix<double,Eigen::Dynamic,Eigen::Dynamic> DirectionX(2,1);
        DirectionX.SetValue(0,0,1.0);

        NuTo::FullMatrix<double,Eigen::Dynamic,Eigen::Dynamic> DirectionY(2,1);
        DirectionY.SetValue(1,0,1.0);

        NuTo::FullMatrix<double,Eigen::Dynamic,Eigen::Dynamic> DirectionXY(2,1);
        DirectionXY.SetValue(0,0,1.0);
        DirectionXY.SetValue(1,0,1.0);






        myStructure.LoadCreateNodeGroupForce(0,GRPLeftBottomSupport,DirectionX, 1);

                                  myStructure.ConstraintLinearSetDisplacementNodeGroup(GRPNodes_Left ,DirectionX,0);
        int ConstraintRightDisp = myStructure.ConstraintLinearSetDisplacementNodeGroup(GRPNodes_Right,DirectionX,0);





        // %%%%%%%%%%%%%
        // Visualization
        // %%%%%%%%%%%%%

#ifdef ENABLE_VISUALIZE
        //mkdir(VTKFolder.c_str(),0777);
        myStructure.AddVisualizationComponentDisplacements();
        myStructure.AddVisualizationComponentEngineeringStrain();
        myStructure.AddVisualizationComponentEngineeringStress();
        //myStructure.ExportVtkDataFileElements(VTKFolder+"/"+VTKFile,false);
#endif





        // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        // Create time integration scheme
        // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        NuTo::NewmarkDirect myTimeIntegrationScheme(&myStructure);

        myTimeIntegrationScheme.SetDampingCoefficientMass(0.05);
        //5 timesteps to capture the quarter wave (5/quarter wave)
        myTimeIntegrationScheme.SetMaxTimeStep(Period/20.);
        myTimeIntegrationScheme.SetTimeStep(SimulationTime/TimeSteps);
        //myTimeIntegrationScheme.SetMaxTimeStep(10);
        myTimeIntegrationScheme.SetMinTimeStep(0.001*myTimeIntegrationScheme.GetMaxTimeStep());


        myTimeIntegrationScheme.SetTimeDependentConstraint(ConstraintRightDisp, DispRHS);



        //set output during the simulation to false
        myStructure.SetShowTime(false);
        myStructure.SetNumProcessors(8);

        myTimeIntegrationScheme.AddResultTime("Time");
        myTimeIntegrationScheme.AddResultGroupNodeForce("Forces_GroupNodes_Left" ,GRPNodes_Left);
        myTimeIntegrationScheme.AddResultGroupNodeForce("Forces_GroupNodes_Right",GRPNodes_Right);

        //set result directory
        bool deleteResultDirectoryFirst(true);
        myTimeIntegrationScheme.SetResultDirectory(VTKFolder,deleteResultDirectoryFirst);

        //solve (perform Newton raphson iteration
        myTimeIntegrationScheme.Solve(SimulationTime);





        myStructure.Info();

    }
    catch (NuTo::Exception e)
    {
        std::cout << e.what() << std::endl;
        return 1;
    }

    std::cout << "Test successful" << std::endl;
    return 0;
}
