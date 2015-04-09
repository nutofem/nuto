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
#include <nuto/mechanics/nodes/NodeCoordinatesDof.h>
#include <nuto/mechanics/structures/unstructured/Structure.h>
#include "nuto/mechanics/timeIntegration/NewmarkDirect.h"

#include <nuto/metamodel/PolynomialLeastSquaresFitting.h>

int main()
{
    try
    {



        // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        // Declaration of neccessary variables
        // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



        bool            EnableSorptionHysteresis            = false;
        bool            EnableModiefiedTangentialStiffness  = false;

        unsigned int    NElements   = 16;                               // Number of elements
        unsigned int    NNodes      = NElements+1;                      // Number of nodes
        double          L           = 0.16;                             // Length of the specimen
        double          Area        = 0.04*0.04;

        double          delta_t     = 1.0/1.0 *     1.0 * 24.0 * 60.0 * 60.0;
        double          t           = 0.0;
        double          t_final     = 1.0/1.0 *   293.0 * 24.0 * 60.0 * 60.0;


        // initial node values
        double          InitialRelativeHumidity         =    0.95;
        double          InitialWaterPhaseFraction       =    0.03;


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
        double          BC_RelativeHumidity             =    0.45;
        double          BC_WaterVolumeFraction;       //=    calculated later
        double          BC_Surface_Moisture_Transfer_RH =    1.0e-10 * delta_t;
        double          BC_Surface_Moisture_Transfer_WVF=    1.0e-7 * delta_t;
        bool            SorptionHistoryDesorption       =    true;

        // sorption hysteresis
        double          Ka                              =    0.26       ;
        double          Kd                              =    0.56       ;
        NuTo::PolynomialLeastSquaresFitting AdsorptionFit;
        NuTo::PolynomialLeastSquaresFitting DesorptionFit;


        // max residual
        double          MaxResidual                     =    1.0e-12 * Area;

        // time measurement
        bool            measureTime                     = true;
        timeval         time_begin, time_end;

        // other
        bool            showPrgress                     = false;
        boost::progress_display ProgressBar(int(t_final/delta_t));

        bool            UseVisualization                = true;
        std::string     VTKFile                         = "ResultsMoistureTransport1D.vtk";
        std::string     VTKFolder                       = "Result";



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


        NuTo::Structure MTStructure1D(1);

        // disable output of calculation times
        MTStructure1D.SetShowTime(false);



        // %%%%%%%%%%%%%%
        // Create Section
        // %%%%%%%%%%%%%%


        int Section1 = MTStructure1D.SectionCreate("Truss");
        MTStructure1D.SectionSetArea(Section1, Area);
        MTStructure1D.SectionSetDOF(Section1,"WaterPhaseFraction RelativeHumidity");



        // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        // Create and set Constitutive Law
        // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


        int ConstLaw = MTStructure1D.ConstitutiveLawCreate("MoistureTransport");

        // set variables
        MTStructure1D.ConstitutiveLawSetAdsorptionCoefficients                                  (ConstLaw,AdsorptionFit.GetPolynomialCoefficients());               // set adsorption coefficients
        MTStructure1D.ConstitutiveLawSetBoundarySurfaceRelativeHumidityTransportCoefficient     (ConstLaw,BC_Surface_Moisture_Transfer_RH);                         // set boundary surface relative humidity transport coefficient
        MTStructure1D.ConstitutiveLawSetBoundarySurfaceWaterVolumeFractionTransportCoefficient  (ConstLaw,BC_Surface_Moisture_Transfer_WVF);                        // set boundary surface water volume fraction transport coefficient
        MTStructure1D.ConstitutiveLawSetDesorptionCoefficients                                  (ConstLaw,DesorptionFit.GetPolynomialCoefficients());               // set desorption coefficients
        MTStructure1D.ConstitutiveLawSetEnableModifiedTangentialStiffness                       (ConstLaw,EnableModiefiedTangentialStiffness);                      // sets whether modified tangential stiffness should be used or not
        MTStructure1D.ConstitutiveLawSetEnableSorptionHysteresis                                (ConstLaw,EnableSorptionHysteresis);                                // sets whether sorption hysteresis should be used or not
        MTStructure1D.ConstitutiveLawSetKa                                                      (ConstLaw,Ka);                                                      // set Ka
        MTStructure1D.ConstitutiveLawSetKd                                                      (ConstLaw,Kd);                                                      // set Kd
        MTStructure1D.ConstitutiveLawSetMassExchangeRate                                        (ConstLaw,MassExchangeRate);                                        // set mass exchange rate
        MTStructure1D.ConstitutiveLawSetPorosity                                                (ConstLaw,Porosity);                                                // set porosity
        MTStructure1D.ConstitutiveLawSetVaporPhaseDiffusionCoefficient                          (ConstLaw,VaporPhaseDiffusionCoefficient);                          // set vapor phase diffusion coefficient
        MTStructure1D.ConstitutiveLawSetVaporPhaseDiffusionExponent                             (ConstLaw,VaporPhaseDiffusionExponent);                             // set vapor phase diffusion exponent
        MTStructure1D.ConstitutiveLawSetVaporPhaseSaturationDensity                             (ConstLaw,VaporPhaseSaturationDensity);                             // set vapor phase saturation density
        MTStructure1D.ConstitutiveLawSetWaterPhaseDensity                                       (ConstLaw,WaterPhaseDensity);                                       // set water phase density
        MTStructure1D.ConstitutiveLawSetWaterPhaseDiffusionCoefficient                          (ConstLaw,WaterPhaseDiffusionCoefficient);                          // set water phase diffusion coefficient
        MTStructure1D.ConstitutiveLawSetWaterPhaseDiffusionExponent                             (ConstLaw,WaterPhaseDiffusionExponent);                             // set water phase diffusion exponent

        // Calculate equilibrium water volume fraction
        InitialWaterPhaseFraction   = MTStructure1D.ConstitutiveLawGetEquilibriumWaterVolumeFraction(ConstLaw,InitialRelativeHumidity,MTStructure1D.ConstitutiveLawGetDesorptionCoefficients(ConstLaw));
        BC_WaterVolumeFraction      = MTStructure1D.ConstitutiveLawGetEquilibriumWaterVolumeFraction(ConstLaw,BC_RelativeHumidity,MTStructure1D.ConstitutiveLawGetDesorptionCoefficients(ConstLaw));

        // %%%%%%%%%%%%
        // Create Nodes
        // %%%%%%%%%%%%


        int nodeNum(0);
        for (unsigned int i=0; i<NNodes; i++)
        {
            NuTo::FullVector<double,Eigen::Dynamic> NodeCoordinates(1);
            NodeCoordinates(0)= i*(L/(NNodes-1));
            MTStructure1D.NodeCreate(nodeNum,std::string("WATERPHASEFRACTION RELATIVEHUMIDITY"),NodeCoordinates,2);
            nodeNum++;
        }

        for (unsigned int i=0; i<NNodes; i++)
        {
            MTStructure1D.NodeGetNodePtr(i)->SetRelativeHumidity(0,InitialRelativeHumidity);
            MTStructure1D.NodeGetNodePtr(i)->SetWaterPhaseFraction(0,InitialWaterPhaseFraction);
        }



        // %%%%%%%%%%%%%%%
        // Create Elements
        // %%%%%%%%%%%%%%%


        NuTo::FullVector<int, Eigen::Dynamic> ElementNodeNumbers(2);

        for (unsigned int i=0; i<NElements; i++)
        {
            // set node numbers
            ElementNodeNumbers(0) = i;
            ElementNodeNumbers(1) = i+1;

            // create elements
            MTStructure1D.ElementCreate(i,"Truss1D2N",ElementNodeNumbers,"CONSTITUTIVELAWIP","STATICDATA");

            // set element section
            MTStructure1D.ElementSetSection(i,Section1);

            // set element constitutive law
            MTStructure1D.ElementSetConstitutiveLaw(i,ConstLaw);
        }



        // %%%%%%%%%%%%%%%%%%%%%%%%
        // Create Boundary Elements
        // %%%%%%%%%%%%%%%%%%%%%%%%


        // Add Nodes to boundary
        int nodeGroupBoundary = MTStructure1D.GroupCreate("NODES");
        MTStructure1D.GroupAddNodeCoordinateRange(nodeGroupBoundary, 0,         -1.e-6,         1.e-6);
        MTStructure1D.GroupAddNodeCoordinateRange(nodeGroupBoundary, 0,         L-1.e-6,        L + 1.e-6);

        // Group all elements with boundary nodes
        int elemGroupBoundary = MTStructure1D.GroupCreate("ELEMENTS");
        MTStructure1D.GroupAddElementsFromNodes(elemGroupBoundary, nodeGroupBoundary, false);

        // Create boundary elements
        ElementNodeNumbers.Resize(1);
        ElementNodeNumbers(0) = 0;
        MTStructure1D.BoundaryElementsCreate("BoundaryMoistureTransport1D", elemGroupBoundary, nodeGroupBoundary, NuTo::BoundaryCondition::NEUMANN_HOMOGENEOUS ,"ConstitutiveLawIp","StaticData");

        // Set Boundary Conditions
        NuTo::FullVector<int,Eigen::Dynamic> BoundaryElementIDs;
        int BoundaryElementGroup = MTStructure1D.GroupCreate("ELEMENTS");
        MTStructure1D.GroupAddElementFromType(BoundaryElementGroup,"BoundaryMoistureTransport1D");
        MTStructure1D.ElementGroupGetMembers(BoundaryElementGroup,BoundaryElementIDs);

        for (int i=0; i<BoundaryElementIDs.GetNumRows(); i++)
        {
            MTStructure1D.ElementGetElementPtr(BoundaryElementIDs(i))->SetBoundaryRelativeHumidity(BC_RelativeHumidity);
            MTStructure1D.ElementGetElementPtr(BoundaryElementIDs(i))->SetBoundaryWaterVolumeFraction(BC_WaterVolumeFraction);
        }



        // %%%%%%%%%%%%%%%
        // Set Static Data
        // %%%%%%%%%%%%%%%


        // Loop over all integration points
        for (int i=0; i<MTStructure1D.GetNumElements(); i++)
        {
            for (int theIP=0; theIP< MTStructure1D.ElementGetElementPtr(i)->GetNumIntegrationPoints(); theIP++)
            {
                NuTo::ConstitutiveStaticDataMoistureTransport *StaticData = MTStructure1D.ElementGetElementPtr(i)->GetStaticData(theIP)->AsMoistureTransport();
                StaticData->SetLastSorptionCoeff(MTStructure1D.ConstitutiveLawGetDesorptionCoefficients(ConstLaw));
                StaticData->SetActualSorptionCoeff(MTStructure1D.ConstitutiveLawGetDesorptionCoefficients(ConstLaw));
                StaticData->SetLastRelHumValue(InitialRelativeHumidity);
                StaticData ->SetSorptionHistoryDesorption(SorptionHistoryDesorption);
            }
        }



        // %%%%%%%%%%%%%%%%%%%%%
        // Multi processor Setup
        // %%%%%%%%%%%%%%%%%%%%%


#ifdef _OPENMP
        MTStructure1D.SetNumProcessors(4);
#else
        MTStructure1D.SetNumProcessors(1);
#endif
        MTStructure1D.CalculateMaximumIndependentSets();
        MTStructure1D.NodeBuildGlobalDofs();



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



        // %%%%%%%%%%%%%%%%%%%%%%
        // Test Getter and Setter
        // %%%%%%%%%%%%%%%%%%%%%%


        // Sorption Curves
        if (AdsorptionFit.GetDegree() != 3)
            throw NuTo::Exception("[Testfile: MoistureTransport1D.cpp] Getter/Setter error: AdsorptionFit.GetDegree() != 3");
        if (DesorptionFit.GetDegree() != 3)
            throw NuTo::Exception("[Testfile: MoistureTransport1D.cpp] Getter/Setter error: DesorptionFit.GetDegree() != 3");

        // Section
        if (MTStructure1D.SectionGetArea(Section1)!=Area)
            throw NuTo::Exception("[Testfile: MoistureTransport1D.cpp] Getter/Setter error: MTStructure1D.SectionGetArea(Section1)!=Area");

        // Constitutive Law
        if (MTStructure1D.ConstitutiveLawGetBoundarySurfaceRelativeHumidityTransportCoefficient(ConstLaw)!=BC_Surface_Moisture_Transfer_RH)
            throw NuTo::Exception("[Testfile: MoistureTransport1D.cpp] Getter/Setter error: MTStructure1D.ConstitutiveLawGetBoundarySurfaceRelativeHumidityTransportCoefficient(ConstLaw)!=BC_Surface_Moisture_Transfer_RH");
        if (MTStructure1D.ConstitutiveLawGetBoundarySurfaceWaterVolumeFractionTransportCoefficient(ConstLaw)!=BC_Surface_Moisture_Transfer_WVF)
            throw NuTo::Exception("[Testfile: MoistureTransport1D.cpp] Getter/Setter error: MTStructure1D.ConstitutiveLawGetBoundarySurfaceWaterVolumeFractionTransportCoefficient(ConstLaw)!=BC_Surface_Moisture_Transfer_WVF");
        if (MTStructure1D.ConstitutiveLawGetEnableModifiedTangentialStiffness(ConstLaw)!=EnableModiefiedTangentialStiffness)
            throw NuTo::Exception("[Testfile: MoistureTransport1D.cpp] Getter/Setter error: MTStructure1D.ConstitutiveLawGetEnableModifiedTangentialStiffness(ConstLaw)!=EnableModiefiedTangentialStiffness");
        if (MTStructure1D.ConstitutiveLawGetEnableSorptionHysteresis(ConstLaw)!=EnableSorptionHysteresis)
            throw NuTo::Exception("[Testfile: MoistureTransport1D.cpp] Getter/Setter error: MTStructure1D.ConstitutiveLawGetEnableSorptionHysteresis(ConstLaw)!=EnableSorptionHysteresis");
        if (MTStructure1D.ConstitutiveLawGetKa(ConstLaw)!=Ka)
            throw NuTo::Exception("[Testfile: MoistureTransport1D.cpp] Getter/Setter error: MTStructure1D.ConstitutiveLawGetKa(ConstLaw)!=Ka");
        if (MTStructure1D.ConstitutiveLawGetKd(ConstLaw)!=Kd)
            throw NuTo::Exception("[Testfile: MoistureTransport1D.cpp] Getter/Setter error: MTStructure1D.ConstitutiveLawGetKd(ConstLaw)!=Kd");
        if (MTStructure1D.ConstitutiveLawGetMassExchangeRate(ConstLaw)!=MassExchangeRate)
            throw NuTo::Exception("[Testfile: MoistureTransport1D.cpp] Getter/Setter error: MTStructure1D.ConstitutiveLawGetMassExchangeRate(ConstLaw)!=MassExchangeRate");
        if (MTStructure1D.ConstitutiveLawGetPorosity(ConstLaw)!=Porosity)
            throw NuTo::Exception("[Testfile: MoistureTransport1D.cpp] Getter/Setter error: MTStructure1D.ConstitutiveLawGetPorosity(ConstLaw)!=Porosity");
        if (MTStructure1D.ConstitutiveLawGetVaporPhaseDiffusionCoefficient(ConstLaw)!=VaporPhaseDiffusionCoefficient)
            throw NuTo::Exception("[Testfile: MoistureTransport1D.cpp] Getter/Setter error: MTStructure1D.ConstitutiveLawGetVaporPhaseDiffusionCoefficient(ConstLaw)!=VaporPhaseDiffusionCoefficient");
        if (MTStructure1D.ConstitutiveLawGetVaporPhaseDiffusionExponent(ConstLaw)!=VaporPhaseDiffusionExponent)
            throw NuTo::Exception("[Testfile: MoistureTransport1D.cpp] Getter/Setter error: MTStructure1D.ConstitutiveLawGetVaporPhaseDiffusionExponent(ConstLaw)!=VaporPhaseDiffusionExponent");
        if (MTStructure1D.ConstitutiveLawGetVaporPhaseSaturationDensity(ConstLaw)!=VaporPhaseSaturationDensity)
            throw NuTo::Exception("[Testfile: MoistureTransport1D.cpp] Getter/Setter error: MTStructure1D.ConstitutiveLawGetVaporPhaseSaturationDensity(ConstLaw)!=VaporPhaseSaturationDensity");
        if (MTStructure1D.ConstitutiveLawGetWaterPhaseDensity(ConstLaw)!=WaterPhaseDensity)
            throw NuTo::Exception("[Testfile: MoistureTransport1D.cpp] Getter/Setter error: MTStructure1D.ConstitutiveLawGetWaterPhaseDensity(ConstLaw)!=WaterPhaseDensity");
        if (MTStructure1D.ConstitutiveLawGetWaterPhaseDiffusionCoefficient(ConstLaw)!=WaterPhaseDiffusionCoefficient)
            throw NuTo::Exception("[Testfile: MoistureTransport1D.cpp] Getter/Setter error: MTStructure1D.ConstitutiveLawGetWaterPhaseDiffusionCoefficient(ConstLaw)!=WaterPhaseDiffusionCoefficient");
        if (MTStructure1D.ConstitutiveLawGetWaterPhaseDiffusionExponent(ConstLaw)!=WaterPhaseDiffusionExponent)
            throw NuTo::Exception("[Testfile: MoistureTransport1D.cpp] Getter/Setter error: MTStructure1D.ConstitutiveLawGetWaterPhaseDiffusionExponent(ConstLaw)!=WaterPhaseDiffusionExponent");

        // Nodes
        for (unsigned int i=0; i<NNodes; i++)
        {
            if (MTStructure1D.NodeGetNodePtr(i)->GetRelativeHumidity(0) != InitialRelativeHumidity)
                throw NuTo::Exception(std::string("[Testfile: MoistureTransport1D.cpp] Getter/Setter error: MTStructure1D.NodeGetNodePtr("+std::to_string(i)+")->GetRelativeHumidity(0) != InitialRelativeHumidity"));
            if (MTStructure1D.NodeGetNodePtr(i)->GetWaterPhaseFraction(0) != InitialWaterPhaseFraction)
                throw NuTo::Exception(std::string("[Testfile: MoistureTransport1D.cpp] Getter/Setter error: MTStructure1D.NodeGetNodePtr("+std::to_string(i)+")->GetRelativeHumidity(0) != InitialRelativeHumidity"));
        }



        // %%%%%%%%%%%%%%%%%%%%%%%
        // Manual time integration
        // %%%%%%%%%%%%%%%%%%%%%%%


        if(measureTime)
        {
            gettimeofday(&time_begin, NULL);
        }

        NuTo::SparseMatrixCSRVector2General<double> Hessian, Hessian0, Hessian1;
        Hessian.Resize(NNodes*2,NNodes*2);
        Hessian0.Resize(NNodes*2,NNodes*2);
        Hessian1.Resize(NNodes*2,NNodes*2);

        NuTo::SparseMatrixCSRGeneral<double> BufferMat;

        NuTo::FullVector<double,Eigen::Dynamic> dis, vel, acc;
        dis.resize(NNodes*2);
        vel.resize(NNodes*2);
        acc.resize(NNodes*2);

        NuTo::FullVector<double,Eigen::Dynamic> dis_last, vel_last, acc_last, delta_dis;
        dis_last.resize(NNodes*2);
        vel_last.resize(NNodes*2);
        acc_last.resize(NNodes*2);
        //delta_dis.resize(NNodes*2);

        NuTo::FullVector<double,Eigen::Dynamic> Buffer_Vec;

        NuTo::FullVector<double,Eigen::Dynamic> res, resComp;
        res.resize(NNodes*2);
        resComp.resize(NNodes*2);

        int NIterations = 0;

        // Sigma-weighted scheme
        // ---------------------

        double sigma = 0.5;


        while (t < t_final)
        {
            NIterations = 0;


            t += delta_t;

            if (!showPrgress)
            {
                std::cout << "Actual timestep: " << t << std::endl;
            }

            MTStructure1D.ElementTotalUpdateStaticData();

            MTStructure1D.NodeExtractDofValues(0,dis_last,Buffer_Vec);
            MTStructure1D.NodeExtractDofValues(1,vel_last,Buffer_Vec);

            // prediction
            //-----------

            vel = vel_last;
            dis = dis_last + (1-sigma) * delta_t * vel_last + sigma * delta_t * vel;

            MTStructure1D.NodeMergeActiveDofValues(0,dis);

            // Calculate residual
            // ------------------
            Buffer_Vec.Resize(0);
            MTStructure1D.BuildGlobalGradientInternalPotentialVector(res);

            do
            {
                NIterations++;

                // calculate correction
                // --------------------

                MTStructure1D.BuildGlobalCoefficientMatrix0(Hessian0, Buffer_Vec);
                MTStructure1D.BuildGlobalCoefficientMatrix1(Hessian1, Buffer_Vec);
                Hessian = Hessian0;
                Hessian.AddScal(Hessian1,1.0/(delta_t * sigma));
                res = -res;

                NuTo::SparseMatrixCSRGeneral<double> HessianForSolver(Hessian);
                HessianForSolver.SetOneBasedIndexing();

                Solver.Solve(HessianForSolver,res,delta_dis);

                // apply correction
                // ----------------

                dis +=  delta_dis;
                vel += 1.0 / (delta_t * sigma) * delta_dis;


                MTStructure1D.NodeMergeActiveDofValues(0,dis);
                MTStructure1D.NodeMergeActiveDofValues(1,vel);


                // Calculate residual
                // ------------------
                MTStructure1D.BuildGlobalGradientInternalPotentialVector(res);

            }
            while(std::abs(res.Min())>MaxResidual || std::abs(res.Max())>MaxResidual);


            if(showPrgress)
            {
                ++ProgressBar;
            }
            else
            {
                std::cout << "Caonvergence after " << NIterations << " iterations "<< std::endl;
                std::cout << "Final residual: " << ((std::abs(res.Min()) > std::abs(res.Max())) ? std::abs(res.Min()) : std::abs(res.Max())) << std::endl << std::endl;
            }
        }

        if(measureTime)
        {
            gettimeofday(&time_end, NULL);
        }



        // %%%%%%%%%%%%%%%%%
        // Visualize Results
        // %%%%%%%%%%%%%%%%%


        if(UseVisualization)
        {
            mkdir(VTKFolder.c_str(),0777);
            MTStructure1D.AddVisualizationComponentRelativeHumidity();
            MTStructure1D.AddVisualizationComponentWaterVolumeFraction();
            MTStructure1D.ExportVtkDataFileElements(VTKFolder+"/"+VTKFile,false);
        }

        // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        // Compare Results with paper from Johannesson and Nyman(2010)
        // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        if(NElements==16)
        {
            NuTo::FullVector<double,Eigen::Dynamic> PaperValues(NNodes);
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

            NuTo::FullVector<double,Eigen::Dynamic> WPF;
            WPF.resize(NNodes);
            for (unsigned int i=0; i<NNodes; i++)
            {
                WPF(i)  = MTStructure1D.NodeGetNodePtr(i)->GetWaterPhaseFraction();
            }

            NuTo::FullVector<double,Eigen::Dynamic> Diff = WPF - PaperValues;
            if(std::abs(Diff.Max()) > 0.01 || std::abs(Diff.Min()) > 0.01)
            {
                throw NuTo::Exception("[Testfile: MoistureTransport1D.cpp]: Results differ to much from given values!");
            }

        }
        else
        {
            throw NuTo::Exception("[Testfile: MoistureTransport1D.cpp]: This testfile needs 16 elements to compare the results with the values taken from the paper of Johannesson and Nyman(2010)");
        }

        if(measureTime)
        {
            std::cout << "elapsed time : " << (time_end.tv_sec - time_begin.tv_sec) + (time_end.tv_usec - time_begin.tv_usec)/1000000.0<< " seconds" << std::endl;
        }
    }    
    catch(NuTo::Exception e)
    {
        std::cout << e.ErrorMessage() << std::endl;
        return 1;
    }
}
