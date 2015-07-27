
#include <fstream>
#include <sys/stat.h>
#include <sys/time.h>

#if defined HAVE_PARDISO
    #include <nuto/math/SparseDirectSolverPardiso.h>
#elif defined HAVE_MUMPS
    #include <nuto/math/SparseDirectSolverMUMPS.h>
#else
    std::cout << "Solver not available - can't solve system of equations " << std::endl;
#endif

#include <boost-1_55/boost/progress.hpp>

#include <nuto/math/SparseMatrixCSRGeneral.h>

#include <nuto/mechanics/constitutive/moistureTransport/ConstitutiveStaticDataMoistureTransport.h>
#include <nuto/mechanics/nodes/NodeDof.h>
#include <nuto/mechanics/nodes/NodeCoordinates.h>
#include <nuto/mechanics/structures/unstructured/Structure.h>
#include "nuto/mechanics/timeIntegration/CrankNicolsonEvaluate.h"

#include <nuto/metamodel/PolynomialLeastSquaresFitting.h>

#include <nuto/mechanics/structures/StructureOutputFullVectorDouble.h>
#include <nuto/mechanics/structures/StructureOutputSparseMatrix.h>

int main()
{
    try
    {



        // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        // Declaration of neccessary variables
        // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



        bool            EnableSorptionHysteresis            = false;
        bool            EnableModiefiedTangentialStiffness  = false;

        //unsigned int    NElements   = 16;                               // Number of elements
        unsigned int    NNodes;                                         // Number of nodes


        double Height           = 0.016;
        double Length           = 0.16;

        unsigned int NumNodesX = 10;
        unsigned int NumNodesY = 2;

        //double ElementHeight    = static_cast<double>(Height)/static_cast<double>(NumNodesY-1);
        //double ElementLength    = static_cast<double>(Length)/static_cast<double>(NumNodesX-1);

        double          Thickness    = 0.04;

        double          delta_t     = 1.0/1.0 *     1.0 * 24.0 * 60.0 * 60.0;
        double          t_final     = 1.0/1.0 *    10.0 * 24.0 * 60.0 * 60.0;


        // initial node values
        double          InitialRelativeHumidity         =    0.95;
        double          InitialWaterVolumeFraction      =    0.03;


        // constitutive law values
        double          MassExchangeRate                =    3.42e-7    ;
        double          Porosity                        =    0.25      ;
        double          VaporPhaseDiffusionCoefficient  =    3.9e-10     ;
        double          VaporPhaseDiffusionExponent     =    1.0        ;
        double          VaporPhaseSaturationDensity     =    0.0173     ;
        double          WaterPhaseDensity               =  999.97       ;
        double          WaterPhaseDiffusionCoefficient  =    1.17e-7    ;
        double          WaterPhaseDiffusionExponent     =    2.0        ;

        // Boundary Condition Values
        //double          BC_RelativeHumidity             =    0.45;
        //double          BC_WaterVolumeFraction;       //=    calculated later
        double          BC_Surface_Moisture_Transfer_RH =    1.0e-10 * delta_t;
        double          BC_Surface_Moisture_Transfer_WVF=    1.0e-7 * delta_t;
        bool            SorptionHistoryDesorption       =    true;

        // sorption hysteresis
        double          Ka                              =    0.26       ;
        double          Kd                              =    0.56       ;
        NuTo::PolynomialLeastSquaresFitting AdsorptionFit;
        NuTo::PolynomialLeastSquaresFitting DesorptionFit;


        // max residual
        double          MaxResidual                     =    1.0e-12;





        // %%%%%%%%%%%%%%%%%%%
        // Fit Sorption Curves
        // %%%%%%%%%%%%%%%%%%%


        NuTo::FullVector <double,Eigen::Dynamic> x_Values_Ad({0.1, 0.2 ,0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9});
        NuTo::FullVector <double,Eigen::Dynamic> y_Values_Ad({0.017, 0.03, 0.04, 0.048, 0.056, 0.066, 0.077, 0.092,0.114});

        NuTo::FullVector <double,Eigen::Dynamic> x_Values_De({0.1, 0.2 ,0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9});
        NuTo::FullVector <double,Eigen::Dynamic> y_Values_De({0.022, 0.039, 0.052, 0.062, 0.072, 0.083, 0.097, 0.118,0.145});

        AdsorptionFit.SetSupportPoints(1,1,x_Values_Ad.Trans(),y_Values_Ad.Trans());
        AdsorptionFit.SetDegree(3);
        AdsorptionFit.AddBoundaryCondition(0.0,0.0);
        AdsorptionFit.AddBoundaryCondition(1.0,0.141);
        AdsorptionFit.BuildDerived();

        DesorptionFit.SetSupportPoints(1,1,x_Values_De.Trans(),y_Values_De.Trans());
        DesorptionFit.SetDegree(3);
        DesorptionFit.AddBoundaryCondition(0.0,0.0);
        DesorptionFit.AddBoundaryCondition(1.0,0.182);
        DesorptionFit.BuildDerived();



        // %%%%%%%%%%%%%%%%
        // Create Structure
        // %%%%%%%%%%%%%%%%


        NuTo::Structure myStructure(2);
        myStructure.SetNumTimeDerivatives(2);

        // disable output of calculation times
        myStructure.SetShowTime(false);
        myStructure.UseMaximumIndependentSets(true);





        // %%%%%%%%%%%%%%
        // Create Section
        // %%%%%%%%%%%%%%


        int mySection = myStructure.SectionCreate("Plane_Stress");
        myStructure.SectionSetThickness(mySection,Thickness);
        myStructure.ElementTotalSetSection(mySection);



        // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        // Create and set Constitutive Law
        // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


        int ConstLaw = myStructure.ConstitutiveLawCreate("MoistureTransport");

        // set variables
        myStructure.ConstitutiveLawSetVariableBool    (ConstLaw,"ENABLE_MODIFIED_TANGENTIAL_STIFFNESS",EnableModiefiedTangentialStiffness);      // sets whether modified tangential stiffness should be used or not
        myStructure.ConstitutiveLawSetVariableBool    (ConstLaw,"enable_sorption_hysteresis",EnableSorptionHysteresis);                          // sets whether sorption hysteresis should be used or not


        myStructure.ConstitutiveLawSetVariableDouble  (ConstLaw,"boundary_TRANSPORT_CONSTANT_GAS_PHASE",BC_Surface_Moisture_Transfer_RH);        // set water phase density
        myStructure.ConstitutiveLawSetVariableDouble  (ConstLaw,"BOUNDARY_TRANSPORT_CONSTANT_WATER_PHASE",BC_Surface_Moisture_Transfer_WVF);     // set water phase density
        myStructure.ConstitutiveLawSetVariableDouble  (ConstLaw,"DENSITY_WATER_PHASE",WaterPhaseDensity);                                        // set water phase density
        myStructure.ConstitutiveLawSetVariableDouble  (ConstLaw,"DIFFUSION_CONSTANT_GAS_PHASE",VaporPhaseDiffusionCoefficient);                  // set vapor phase diffusion coefficient
        myStructure.ConstitutiveLawSetVariableDouble  (ConstLaw,"DIFFUSION_CONSTANT_WATER_PHASE",WaterPhaseDiffusionCoefficient);                // set water phase diffusion coefficient
        myStructure.ConstitutiveLawSetVariableDouble  (ConstLaw,"DIFFUSION_EXPONENT_GAS_PHASE",VaporPhaseDiffusionExponent);                     // set vapor phase diffusion exponent
        myStructure.ConstitutiveLawSetVariableDouble  (ConstLaw,"DIFFUSION_EXPONENT_WATER_PHASE",WaterPhaseDiffusionExponent);                   // set water phase diffusion exponent
        myStructure.ConstitutiveLawSetVariableDouble  (ConstLaw,"GRADIENT_CORRECTION_ADSORPTION_DESORPTION",Kd);                                 // set gradient correction when changing from adsorption to desorption
        myStructure.ConstitutiveLawSetVariableDouble  (ConstLaw,"GRADIENT_CORRECTION_DESORPTION_ADSORPTION",Ka);                                 // set gradient correction when changing from desorption to adsorption
        myStructure.ConstitutiveLawSetVariableDouble  (ConstLaw,"MASS_EXCHANGE_RATE",MassExchangeRate);                                          // set mass exchange rate
        myStructure.ConstitutiveLawSetVariableDouble  (ConstLaw,"POROSITY",Porosity);                                                            // set porosity
        myStructure.ConstitutiveLawSetVariableDouble  (ConstLaw,"SATURATION_DENSITY_GAS_PHASE",VaporPhaseSaturationDensity);                     // set vapor phase saturation density

        myStructure.ConstitutiveLawSetVariableFullVectorDouble    (ConstLaw,"polynomial_COEFFICIENTS_ADSORPTION",AdsorptionFit.GetPolynomialCoefficients());               // set adsorption coefficients
        myStructure.ConstitutiveLawSetVariableFullVectorDouble    (ConstLaw,"POLYNOMIAL_COEFFICIENTS_DESORPTION",DesorptionFit.GetPolynomialCoefficients());               // set desorption coefficients


        // Calculate equilibrium water volume fraction
        InitialWaterVolumeFraction   = myStructure.ConstitutiveLawGetEquilibriumWaterVolumeFraction(ConstLaw,InitialRelativeHumidity,myStructure.ConstitutiveLawGetVariableFullVectorDouble(ConstLaw,NuTo::Constitutive::eConstitutiveVariable::POLYNOMIAL_COEFFICIENTS_DESORPTION));
        //BC_WaterVolumeFraction       = myStructure.ConstitutiveLawGetEquilibriumWaterVolumeFraction(ConstLaw,BC_RelativeHumidity,myStructure.ConstitutiveLawGetVariableFullVectorDouble(ConstLaw,NuTo::Constitutive::eConstitutiveVariable::POLYNOMIAL_COEFFICIENTS_DESORPTION));



        // %%%%%%%%%%%%
        // Create Nodes
        // %%%%%%%%%%%%




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
        myStructure.InterpolationTypeAdd(myInterpolationType, NuTo::Node::RELATIVEHUMIDITY, NuTo::Interpolation::eTypeOrder::EQUIDISTANT1);
        myStructure.InterpolationTypeAdd(myInterpolationType, NuTo::Node::WATERVOLUMEFRACTION, NuTo::Interpolation::eTypeOrder::EQUIDISTANT1);






        // %%%%%%%%%%%%%%%
        // Create Elements
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
                int elementNumber = myStructure.ElementCreate(myInterpolationType, Nodes, "CONSTITUTIVELAWIP","STATICDATA");

                // set element section
                myStructure.ElementSetSection(elementNumber,mySection);

                // set element constitutive law
                myStructure.ElementSetConstitutiveLaw(elementNumber,ConstLaw);
            }
        }

        myStructure.ElementTotalConvertToInterpolationType();

//        // Apply nodal start values

        NNodes = myStructure.GetNumNodes();

        for (unsigned int i=0; i<NNodes; i++)
        {
            auto NodeMultiplier= myStructure.NodeGetNodePtr(i)->GetCoordinates()[0]/Length;
            if(myStructure.NodeGetNodePtr(i)->GetNumRelativeHumidity() != 0)
            {
                myStructure.NodeGetNodePtr(i)->SetRelativeHumidity(0,InitialRelativeHumidity*NodeMultiplier);
            }
            if(myStructure.NodeGetNodePtr(i)->GetNumWaterVolumeFraction() != 0)
            {
                myStructure.NodeGetNodePtr(i)->SetWaterVolumeFraction(0,InitialWaterVolumeFraction*NodeMultiplier);
            }
        }

//        // %%%%%%%%%%%%%%%%%%%%%%%%
//        // Create Boundary Elements
//        // %%%%%%%%%%%%%%%%%%%%%%%%


//        // Add Nodes to boundary
//        int nodeGroupBoundary = myStructure.GroupCreate("NODES");
//        myStructure.GroupAddNodeCoordinateRange(nodeGroupBoundary, 0,         -1.e-6,         1.e-6);
//        myStructure.GroupAddNodeCoordinateRange(nodeGroupBoundary, 0,         L-1.e-6,        L + 1.e-6);

//        // Group all elements with boundary nodes
//        int elemGroupBoundary = myStructure.GroupCreate("ELEMENTS");
//        myStructure.GroupAddElementsFromNodes(elemGroupBoundary, nodeGroupBoundary, false);

//        // Create boundary elements
//        ElementNodeNumbers.Resize(1);
//        ElementNodeNumbers(0) = 0;

//        auto BoundaryElementIDs = myStructure.BoundaryElementsCreate(elemGroupBoundary,nodeGroupBoundary);

//        for (int i=0; i<BoundaryElementIDs.GetNumRows(); i++)
//        {
//            myStructure.ElementGetElementPtr(BoundaryElementIDs(i))->SetBoundaryRelativeHumidity(BC_RelativeHumidity);
//            myStructure.ElementGetElementPtr(BoundaryElementIDs(i))->SetBoundaryWaterVolumeFraction(BC_WaterVolumeFraction);
//        }

//        // %%%%%%%%%%%%%%%
//        // Set Static Data
//        // %%%%%%%%%%%%%%%


        // Loop over all integration points
        for (int i=0; i<myStructure.GetNumElements(); i++)
        {
            for (int theIP=0; theIP< myStructure.ElementGetElementPtr(i)->GetNumIntegrationPoints(); theIP++)
            {
                NuTo::ConstitutiveStaticDataMoistureTransport *StaticData = myStructure.ElementGetElementPtr(i)->GetStaticData(theIP)->AsMoistureTransport();
                StaticData->SetLastSorptionCoeff(myStructure.ConstitutiveLawGetVariableFullVectorDouble(ConstLaw,NuTo::Constitutive::eConstitutiveVariable::POLYNOMIAL_COEFFICIENTS_DESORPTION));
                StaticData->SetActualSorptionCoeff(myStructure.ConstitutiveLawGetVariableFullVectorDouble(ConstLaw,NuTo::Constitutive::eConstitutiveVariable::POLYNOMIAL_COEFFICIENTS_DESORPTION));
                StaticData->SetLastRelHumValue(InitialRelativeHumidity);
                StaticData ->SetSorptionHistoryDesorption(SorptionHistoryDesorption);
            }
        }



        // %%%%%%%%%%%%%%%%%%%%%
        // Multi processor Setup
        // %%%%%%%%%%%%%%%%%%%%%


#ifdef _OPENMP
        myStructure.SetNumProcessors(4);
        std::cout << "OpenMP enabled" << std::endl;
#else
        myStructure.SetNumProcessors(1);
#endif
        myStructure.CalculateMaximumIndependentSets();
        myStructure.UseMaximumIndependentSets(true);
        myStructure.NodeBuildGlobalDofs();



        //myStructure.NodeInfo(10);

        // %%%%%%%%%%
        // Set Solver
        // %%%%%%%%%%


#if defined HAVE_PARDISO
        NuTo::SparseDirectSolverPardiso Solver(4);
#elif defined HAVE_MUMPS
        NuTo::SparseDirectSolverMUMPS Solver;
#else
        std::cout << "Solver not available - can't solve system of equations " << std::endl;
#endif
        Solver.SetShowTime(false);


        // %%%%%%%%%%%%%%%%
        // Time Integration
        // %%%%%%%%%%%%%%%%

        myStructure.AddVisualizationComponentRelativeHumidity();
        myStructure.AddVisualizationComponentWaterVolumeFraction();

        NuTo::CrankNicolsonEvaluate TimeIntegrationScheme(&myStructure);

        TimeIntegrationScheme.SetTimeStep(delta_t);
        TimeIntegrationScheme.SetMaxTimeStep(delta_t);
        TimeIntegrationScheme.SetMinTimeStep(delta_t);

        TimeIntegrationScheme.SetPerformLineSearch(false);
        TimeIntegrationScheme.SetCheckEquilibriumOnStart(false);
        TimeIntegrationScheme.SetToleranceForce(MaxResidual);
        TimeIntegrationScheme.SetMaxNumIterations(40);

        //set result directory
        bool deleteResultDirectoryFirst(true);
        TimeIntegrationScheme.SetResultDirectory("./ResultsMoistureTransport1D",deleteResultDirectoryFirst);

        TimeIntegrationScheme.Solve(t_final);






//        // %%%%%%%%%%%%%%%%%
//        // Visualize Results
//        // %%%%%%%%%%%%%%%%%


//        if(UseVisualization)
//        {
//            mkdir(VTKFolder.c_str(),0777);
//            myStructure.ExportVtkDataFileElements(VTKFolder+"/"+VTKFile,false);
//        }

//        // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//        // Compare Results with paper from Johannesson and Nyman(2010)
//        // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

//        // Just in case one uses an interpolation order higher than 1
//        NNodes = NElements + 1;

//        if(NElements==16)
//        {
//            NuTo::FullVector<double,Eigen::Dynamic> PaperValues(NNodes);
//            PaperValues(0)  = 0.06;
//            PaperValues(1)  = 0.097;
//            PaperValues(2)  = 0.116;
//            PaperValues(3)  = 0.129;
//            PaperValues(4)  = 0.138;
//            PaperValues(5)  = 0.146;
//            PaperValues(6)  = 0.148;
//            PaperValues(7)  = 0.151;
//            PaperValues(8)  = 0.152;
//            PaperValues(9)  = PaperValues(7);
//            PaperValues(10) = PaperValues(6);
//            PaperValues(11) = PaperValues(5);
//            PaperValues(12) = PaperValues(4);
//            PaperValues(13) = PaperValues(3);
//            PaperValues(14) = PaperValues(2);
//            PaperValues(15) = PaperValues(1);
//            PaperValues(16) = PaperValues(0);

//            NuTo::FullVector<double,Eigen::Dynamic> WPF;
//            WPF.resize(NNodes);
//            for (unsigned int i=0; i<NNodes; i++)
//            {
//                WPF(i)  = myStructure.NodeGetNodePtr(i)->GetWaterVolumeFraction();
//            }

//            NuTo::FullVector<double,Eigen::Dynamic> Diff = WPF - PaperValues;
//            if(std::abs(Diff.Max()) > 0.01 || std::abs(Diff.Min()) > 0.01)
//            {
//                throw NuTo::Exception("[Testfile: MoistureTransport1D.cpp]: Results differ to much from given values!");
//            }
//            else
//            {
//                std::cout << "Calculation finished!" << std::endl;
//            }

//        }
//        else
//        {
//            throw NuTo::Exception("[Testfile: MoistureTransport1D.cpp]: This testfile needs 16 elements to compare the results with the values taken from the paper of Johannesson and Nyman(2010)");
//        }

//        if(measureTime)
//        {
//            std::cout << "elapsed time : " << (time_end.tv_sec - time_begin.tv_sec) + (time_end.tv_usec - time_begin.tv_usec)/1000000.0<< " seconds" << std::endl;
//        }
    }
    catch(NuTo::Exception e)
    {
        std::cout << e.ErrorMessage() << std::endl;
        return 1;
    }
}
