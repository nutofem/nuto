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

#include <boost/progress.hpp>

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

        unsigned int    NElements   = 16;                               // Number of elements
        unsigned int    NNodes      = NElements+1;                      // Number of nodes
        double          L           = 0.16;                             // Length of the specimen
        double          Area        = 0.04*0.04;

        double          delta_t     = 1.0/1.0 *     1.0 * 24.0 * 60.0 * 60.0;
        double          t_final     = 1.0/1.0 *   293.0 * 24.0 * 60.0 * 60.0;


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
        double          MaxResidual                     =    1.0e-13;

        // time measurement
        bool            measureTime                     = true;
        timeval         time_begin, time_end;



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
        MTStructure1D.SetNumTimeDerivatives(2);

        // disable output of calculation times
        MTStructure1D.SetShowTime(false);



        // %%%%%%%%%%%%%%
        // Create Section
        // %%%%%%%%%%%%%%


        int Section1 = MTStructure1D.SectionCreate("Truss");
        MTStructure1D.SectionSetArea(Section1, Area);



        // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        // Create and set Constitutive Law
        // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


        int ConstLaw = MTStructure1D.ConstitutiveLawCreate("MoistureTransport");

        // set variables
        MTStructure1D.ConstitutiveLawSetParameterBool    (ConstLaw,"ENABLE_MODIFIED_TANGENTIAL_STIFFNESS",EnableModiefiedTangentialStiffness);      // sets whether modified tangential stiffness should be used or not
        MTStructure1D.ConstitutiveLawSetParameterBool    (ConstLaw,"enable_sorption_hysteresis",EnableSorptionHysteresis);                          // sets whether sorption hysteresis should be used or not


        MTStructure1D.ConstitutiveLawSetParameterDouble  (ConstLaw,"boundary_TRANSPORT_CONSTANT_GAS_PHASE",BC_Surface_Moisture_Transfer_RH);        // set water phase density
        MTStructure1D.ConstitutiveLawSetParameterDouble  (ConstLaw,"BOUNDARY_TRANSPORT_CONSTANT_WATER_PHASE",BC_Surface_Moisture_Transfer_WVF);     // set water phase density
        MTStructure1D.ConstitutiveLawSetParameterDouble  (ConstLaw,"DENSITY_WATER_PHASE",WaterPhaseDensity);                                        // set water phase density
        MTStructure1D.ConstitutiveLawSetParameterDouble  (ConstLaw,"DIFFUSION_CONSTANT_GAS_PHASE",VaporPhaseDiffusionCoefficient);                  // set vapor phase diffusion coefficient
        MTStructure1D.ConstitutiveLawSetParameterDouble  (ConstLaw,"DIFFUSION_CONSTANT_WATER_PHASE",WaterPhaseDiffusionCoefficient);                // set water phase diffusion coefficient
        MTStructure1D.ConstitutiveLawSetParameterDouble  (ConstLaw,"DIFFUSION_EXPONENT_GAS_PHASE",VaporPhaseDiffusionExponent);                     // set vapor phase diffusion exponent
        MTStructure1D.ConstitutiveLawSetParameterDouble  (ConstLaw,"DIFFUSION_EXPONENT_WATER_PHASE",WaterPhaseDiffusionExponent);                   // set water phase diffusion exponent
        MTStructure1D.ConstitutiveLawSetParameterDouble  (ConstLaw,"GRADIENT_CORRECTION_ADSORPTION_DESORPTION",Kd);                                 // set gradient correction when changing from adsorption to desorption
        MTStructure1D.ConstitutiveLawSetParameterDouble  (ConstLaw,"GRADIENT_CORRECTION_DESORPTION_ADSORPTION",Ka);                                 // set gradient correction when changing from desorption to adsorption
        MTStructure1D.ConstitutiveLawSetParameterDouble  (ConstLaw,"MASS_EXCHANGE_RATE",MassExchangeRate);                                          // set mass exchange rate
        MTStructure1D.ConstitutiveLawSetParameterDouble  (ConstLaw,"POROSITY",Porosity);                                                            // set porosity
        MTStructure1D.ConstitutiveLawSetParameterDouble  (ConstLaw,"SATURATION_DENSITY_GAS_PHASE",VaporPhaseSaturationDensity);                     // set vapor phase saturation density

        MTStructure1D.ConstitutiveLawSetParameterFullVectorDouble    (ConstLaw,"polynomial_COEFFICIENTS_ADSORPTION",AdsorptionFit.GetPolynomialCoefficients());               // set adsorption coefficients
        MTStructure1D.ConstitutiveLawSetParameterFullVectorDouble    (ConstLaw,"POLYNOMIAL_COEFFICIENTS_DESORPTION",DesorptionFit.GetPolynomialCoefficients());               // set desorption coefficients


        // Calculate equilibrium water volume fraction
        InitialWaterVolumeFraction   = MTStructure1D.ConstitutiveLawGetEquilibriumWaterVolumeFraction(ConstLaw,InitialRelativeHumidity,MTStructure1D.ConstitutiveLawGetParameterFullVectorDouble(ConstLaw,NuTo::Constitutive::eConstitutiveParameter::POLYNOMIAL_COEFFICIENTS_DESORPTION));
        BC_WaterVolumeFraction      = MTStructure1D.ConstitutiveLawGetEquilibriumWaterVolumeFraction(ConstLaw,BC_RelativeHumidity,MTStructure1D.ConstitutiveLawGetParameterFullVectorDouble(ConstLaw,NuTo::Constitutive::eConstitutiveParameter::POLYNOMIAL_COEFFICIENTS_DESORPTION));



        // %%%%%%%%%%%%
        // Create Nodes
        // %%%%%%%%%%%%


        int nodeNum(0);
        for (unsigned int i=0; i<NNodes; i++)
        {
            NuTo::FullVector<double,Eigen::Dynamic> NodeCoordinates(1);
            NodeCoordinates(0)= i*(L/(NNodes-1));
            MTStructure1D.NodeCreate(nodeNum,NodeCoordinates);
            nodeNum++;
        }



        // %%%%%%%%%%%%%%%%%%
        // Interpolation Type
        // %%%%%%%%%%%%%%%%%%

        int Interpol = MTStructure1D.InterpolationTypeCreate("Truss1D");
        MTStructure1D.InterpolationTypeAdd(Interpol, NuTo::Node::COORDINATES, NuTo::Interpolation::eTypeOrder::EQUIDISTANT1);
        MTStructure1D.InterpolationTypeAdd(Interpol, NuTo::Node::RELATIVEHUMIDITY, NuTo::Interpolation::eTypeOrder::EQUIDISTANT1);
        MTStructure1D.InterpolationTypeAdd(Interpol, NuTo::Node::WATERVOLUMEFRACTION, NuTo::Interpolation::eTypeOrder::EQUIDISTANT1);





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
            MTStructure1D.ElementCreate(i,Interpol,ElementNodeNumbers,"CONSTITUTIVELAWIP","STATICDATA");

            // set element section
            MTStructure1D.ElementSetSection(i,Section1);

            // set element constitutive law
            MTStructure1D.ElementSetConstitutiveLaw(i,ConstLaw);
        }

        MTStructure1D.ElementTotalConvertToInterpolationType();

        // Apply nodal start values

        NNodes = MTStructure1D.GetNumNodes();

        for (unsigned int i=0; i<NNodes; i++)
        {
            auto NodeMultiplier= 1.0;//MTStructure1D.NodeGetNodePtr(i)->GetCoordinates()[0]/L;
            if(MTStructure1D.NodeGetNodePtr(i)->GetNumRelativeHumidity() != 0)
            {
                MTStructure1D.NodeGetNodePtr(i)->SetRelativeHumidity(0,InitialRelativeHumidity*NodeMultiplier);
            }
            if(MTStructure1D.NodeGetNodePtr(i)->GetNumWaterVolumeFraction() != 0)
            {
                MTStructure1D.NodeGetNodePtr(i)->SetWaterVolumeFraction(0,InitialWaterVolumeFraction*NodeMultiplier);
            }
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

        int groupBoundaryElements = MTStructure1D.BoundaryElementsCreate(elemGroupBoundary,nodeGroupBoundary);
        auto BoundaryElementIDs = MTStructure1D.GroupGetMemberIds(groupBoundaryElements);
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
                StaticData->SetLastSorptionCoeff(MTStructure1D.ConstitutiveLawGetParameterFullVectorDouble(ConstLaw,NuTo::Constitutive::eConstitutiveParameter::POLYNOMIAL_COEFFICIENTS_DESORPTION));
                StaticData->SetActualSorptionCoeff(MTStructure1D.ConstitutiveLawGetParameterFullVectorDouble(ConstLaw,NuTo::Constitutive::eConstitutiveParameter::POLYNOMIAL_COEFFICIENTS_DESORPTION));
                StaticData->SetLastRelHumValue(InitialRelativeHumidity);
                StaticData ->SetSorptionHistoryDesorption(SorptionHistoryDesorption);
            }
        }



        // %%%%%%%%%%%%%%%%%%%%%
        // Multi processor Setup
        // %%%%%%%%%%%%%%%%%%%%%


#ifdef _OPENMP
        MTStructure1D.SetNumProcessors(4);
        std::cout << "OpenMP enabled" << std::endl;
#else
        MTStructure1D.SetNumProcessors(1);
#endif
        MTStructure1D.CalculateMaximumIndependentSets();
        MTStructure1D.UseMaximumIndependentSets(true);
        MTStructure1D.NodeBuildGlobalDofs();



        //MTStructure1D.NodeInfo(10);

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
        if (MTStructure1D.ConstitutiveLawGetParameterBool(ConstLaw,"enable_MODIFIED_TANGENTIAL_STIFFNESS")!=EnableModiefiedTangentialStiffness)
            throw NuTo::Exception("[Testfile: MoistureTransport1D.cpp] Getter/Setter error: MTStructure1D.ConstitutiveLawGetParameterBool(ConstLaw,\"enable_MODIFIED_TANGENTIAL_STIFFNESS\")!=EnableModiefiedTangentialStiffness");
        if (MTStructure1D.ConstitutiveLawGetParameterBool(ConstLaw,"ENABLE_SORPTION_HYSTERESIS")!=EnableSorptionHysteresis)
            throw NuTo::Exception("[Testfile: MoistureTransport1D.cpp] Getter/Setter error: MTStructure1D.ConstitutiveLawGetParameterBool(ConstLaw,\"ENABLE_SORPTION_HYSTERESIS\")!=EnableSorptionHysteresis");

        if (MTStructure1D.ConstitutiveLawGetParameterDouble(ConstLaw,"boundary_TRANSPORT_CONSTANT_GAS_PHASE")!=BC_Surface_Moisture_Transfer_RH)
            throw NuTo::Exception("[Testfile: MoistureTransport1D.cpp] Getter/Setter error: MTStructure1D.ConstitutiveLawGetParameterDouble(ConstLaw,\"boundary_TRANSPORT_CONSTANT_GAS_PHASE\")!=BC_Surface_Moisture_Transfer_RH");
        if (MTStructure1D.ConstitutiveLawGetParameterDouble(ConstLaw,"BOUNDARY_TRANSPORT_CONSTANT_WATER_PHASE")!=BC_Surface_Moisture_Transfer_WVF)
            throw NuTo::Exception("[Testfile: MoistureTransport1D.cpp] Getter/Setter error: MTStructure1D.ConstitutiveLawGetParameterDouble(ConstLaw,\"BOUNDARY_TRANSPORT_CONSTANT_WATER_PHASE\")!=BC_Surface_Moisture_Transfer_WVF");
        if (MTStructure1D.ConstitutiveLawGetParameterDouble(ConstLaw,"DENSITY_WATER_PHASE")!=WaterPhaseDensity)
            throw NuTo::Exception("[Testfile: MoistureTransport1D.cpp] Getter/Setter error: MTStructure1D.ConstitutiveLawGetParameterDouble(ConstLaw,\"DENSITY_WATER_PHASE\")!=WaterPhaseDensity");
        if (MTStructure1D.ConstitutiveLawGetParameterDouble(ConstLaw,"DIFFUSION_CONSTANT_GAS_PHASE")!=VaporPhaseDiffusionCoefficient)
            throw NuTo::Exception("[Testfile: MoistureTransport1D.cpp] Getter/Setter error: MTStructure1D.ConstitutiveLawGetParameterDouble(ConstLaw,\"DIFFUSION_CONSTANT_GAS_PHASE\")!=VaporPhaseDiffusionCoefficient");
        if (MTStructure1D.ConstitutiveLawGetParameterDouble(ConstLaw,"DIFFUSION_CONSTANT_WATER_PHASE")!=WaterPhaseDiffusionCoefficient)
            throw NuTo::Exception("[Testfile: MoistureTransport1D.cpp] Getter/Setter error: MTStructure1D.ConstitutiveLawGetParameterDouble(ConstLaw,\"DIFFUSION_CONSTANT_WATER_PHASE\")!=WaterPhaseDiffusionCoefficient");
        if (MTStructure1D.ConstitutiveLawGetParameterDouble(ConstLaw,"DIFFUSION_EXPONENT_GAS_PHASE")!=VaporPhaseDiffusionExponent)
            throw NuTo::Exception("[Testfile: MoistureTransport1D.cpp] Getter/Setter error: MTStructure1D.ConstitutiveLawGetParameterDouble(ConstLaw,\"DIFFUSION_EXPONENT_GAS_PHASE\")!=VaporPhaseDiffusionExponent");
        if (MTStructure1D.ConstitutiveLawGetParameterDouble(ConstLaw,"DIFFUSION_EXPONENT_WATER_PHASE")!=WaterPhaseDiffusionExponent)
            throw NuTo::Exception("[Testfile: MoistureTransport1D.cpp] Getter/Setter error: MTStructure1D.ConstitutiveLawGetParameterDouble(ConstLaw,\"DIFFUSION_EXPONENT_WATER_PHASE\")!=WaterPhaseDiffusionExponent");
        if (MTStructure1D.ConstitutiveLawGetParameterDouble(ConstLaw,"GRADIENT_CORRECTION_DESORPTION_ADSORPTION")!=Ka)
            throw NuTo::Exception("[Testfile: MoistureTransport1D.cpp] Getter/Setter error: MTStructure1D.ConstitutiveLawGetParameterDouble(ConstLaw,\"GRADIENT_CORRECTION_DESORPTION_ADSORPTION\")!=Ka");
        if (MTStructure1D.ConstitutiveLawGetParameterDouble(ConstLaw,"GRADIENT_CORRECTION_ADSORPTION_DESORPTION")!=Kd)
            throw NuTo::Exception("[Testfile: MoistureTransport1D.cpp] Getter/Setter error: MTStructure1D.ConstitutiveLawGetParameterDouble(ConstLaw,\"GRADIENT_CORRECTION_ADSORPTION_DESORPTION\")!=Kd");
        if (MTStructure1D.ConstitutiveLawGetParameterDouble(ConstLaw,"MASS_EXCHANGE_RATE")!=MassExchangeRate)
            throw NuTo::Exception("[Testfile: MoistureTransport1D.cpp] Getter/Setter error: MTStructure1D.ConstitutiveLawGetParameterDouble(ConstLaw,\"MASS_EXCHANGE_RATE\")!=MassExchangeRate");
        if (MTStructure1D.ConstitutiveLawGetParameterDouble(ConstLaw,"POROSITY")!=Porosity)
            throw NuTo::Exception("[Testfile: MoistureTransport1D.cpp] Getter/Setter error: MTStructure1D.ConstitutiveLawGetParameterDouble(ConstLaw,\"POROSITY\")!=Porosity");
        if (MTStructure1D.ConstitutiveLawGetParameterDouble(ConstLaw,"SATURATION_DENSITY_GAS_PHASE")!=VaporPhaseSaturationDensity)
            throw NuTo::Exception("[Testfile: MoistureTransport1D.cpp] Getter/Setter error: MTStructure1D.ConstitutiveLawGetParameterDouble(ConstLaw,\"SATURATION_DENSITY_GAS_PHASE\")!=VaporPhaseSaturationDensity");

        // Nodes
        for (unsigned int i=0; i<NNodes; i++)
        {
            if (MTStructure1D.NodeGetNodePtr(i)->GetNumRelativeHumidity() != 0 && MTStructure1D.NodeGetNodePtr(i)->GetRelativeHumidity(0) != InitialRelativeHumidity)
                throw NuTo::Exception(std::string("[Testfile: MoistureTransport1D.cpp] Getter/Setter error: MTStructure1D.NodeGetNodePtr("+std::to_string(i)+")->GetRelativeHumidity(0) != InitialRelativeHumidity"));
            if (MTStructure1D.NodeGetNodePtr(i)->GetNumWaterVolumeFraction() != 0 && MTStructure1D.NodeGetNodePtr(i)->GetWaterVolumeFraction(0) != InitialWaterVolumeFraction)
                throw NuTo::Exception(std::string("[Testfile: MoistureTransport1D.cpp] Getter/Setter error: MTStructure1D.NodeGetNodePtr("+std::to_string(i)+")->GetRelativeHumidity(0) != InitialRelativeHumidity"));
        }

        // Begin Test Structure Evaluate
        // -----------------------------

/*
        // %%%%%%%%%%%%%%%%%%%%%%%
        // Manual time integration
        // %%%%%%%%%%%%%%%%%%%%%%%


        if(measureTime)
        {
            gettimeofday(&time_begin, NULL);
        }

        // Evaluate Version


        NuTo::SparseMatrixCSRVector2General<double> Hessian_0;
        NuTo::SparseMatrixCSRVector2General<double> Hessian_1;

        NuTo::FullVector<double, Eigen::Dynamic> Residual;

        NuTo::StructureOutputFullVectorDouble SOResidual(Residual);
        NuTo::StructureOutputSparseMatrix SOHessian_0(Hessian_0);
        NuTo::StructureOutputSparseMatrix SOHessian_1(Hessian_1);

        std::map<NuTo::StructureEnum::eOutput, NuTo::StructureOutputBase*> StructureOutputResidual;
        std::map<NuTo::StructureEnum::eOutput, NuTo::StructureOutputBase*> StructureOutputHessian;
        std::map<NuTo::StructureEnum::eOutput, NuTo::StructureOutputBase*> StructureOutputResAndHessian;

        StructureOutputResidual[NuTo::StructureEnum::eOutput::INTERNAL_GRADIENT]  = &SOResidual;

        StructureOutputHessian[NuTo::StructureEnum::eOutput::STIFFNESS]  = &SOHessian_0;
        StructureOutputHessian[NuTo::StructureEnum::eOutput::DAMPING]    = &SOHessian_1;

        StructureOutputResAndHessian[NuTo::StructureEnum::eOutput::INTERNAL_GRADIENT]   = &SOResidual;
        StructureOutputResAndHessian[NuTo::StructureEnum::eOutput::STIFFNESS]      = &SOHessian_0;
        StructureOutputResAndHessian[NuTo::StructureEnum::eOutput::DAMPING]        = &SOHessian_1;


        NuTo::SparseMatrixCSRVector2General<double> Hessian;

        //NuTo::SparseMatrixCSRGeneral<double> BufferMat;

        NuTo::FullVector<double,Eigen::Dynamic> dis(NNodes*2), vel(NNodes*2);


        NuTo::FullVector<double,Eigen::Dynamic> dis_last(NNodes*2), vel_last(NNodes*2), delta_dis(NNodes*2);
        //dis_last.resize(NNodes*2);
        //vel_last.resize(NNodes*2);
        //acc_last.resize(NNodes*2);
        //delta_dis.resize(NNodes*2);

        NuTo::FullVector<double,Eigen::Dynamic> Buffer_Vec;


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
            //MTStructure1D.BuildGlobalGradientInternalPotentialVector(Residual);

            //res.Info(12,5,true);

            MTStructure1D.Evaluate(StructureOutputResAndHessian);

            //MTStructure1D.Evaluate(StructureOutputResidual);
            //Residual.GetFullVectorDouble().Info(12,5,true);


            do
            {
                NIterations++;

                // calculate correction
                // --------------------



                //MTStructure1D.BuildGlobalCoefficientMatrix0(Hessian0, Buffer_Vec);
                //MTStructure1D.BuildGlobalCoefficientMatrix1(Hessian1, Buffer_Vec);
                //Hessian = Hessian0;
                //Hessian.AddScal(Hessian1,1.0/(delta_t * sigma));
                if(NIterations > 1)
                {
                    MTStructure1D.Evaluate(StructureOutputHessian);
                }
                Hessian = Hessian_0;
                Hessian.AddScal(Hessian_1,1.0/(delta_t * sigma));

                //res = -Residual.GetFullVectorDouble();

                NuTo::SparseMatrixCSRGeneral<double> HessianForSolver(Hessian);
                HessianForSolver.SetOneBasedIndexing();

                Solver.Solve(HessianForSolver,-Residual,delta_dis);

                // apply correction
                // ----------------

                dis +=  delta_dis;
                vel += 1.0 / (delta_t * sigma) * delta_dis;


                MTStructure1D.NodeMergeActiveDofValues(0,dis);
                MTStructure1D.NodeMergeActiveDofValues(1,vel);


                // Calculate residual
                // ------------------
                //MTStructure1D.BuildGlobalGradientInternalPotentialVector(Residual);
                MTStructure1D.Evaluate(StructureOutputResidual);


            }
            while(std::abs(Residual.Min())>MaxResidual || std::abs(Residual.Max())>MaxResidual);



            if(showPrgress)
            {
                ++ProgressBar;
            }
            else
            {
                std::cout << "Caonvergence after " << NIterations << " iterations "<< std::endl;
                std::cout << "Final residual: " << ((std::abs(Residual.Min()) > std::abs(Residual.Max())) ? std::abs(Residual.Min()) : std::abs(Residual.Max())) << std::endl << std::endl;
            }
        }


        //dis.Info(12,5,true);
*/
#ifdef ENABLE_VISUALIZE
        int visualizationGroup = MTStructure1D.GroupCreate(NuTo::Groups::eGroupId::Elements);
        MTStructure1D.GroupAddElementsTotal(visualizationGroup);

        MTStructure1D.AddVisualizationComponent(visualizationGroup, NuTo::VisualizeBase::RELATIVE_HUMIDITY);
        MTStructure1D.AddVisualizationComponent(visualizationGroup, NuTo::VisualizeBase::WATER_VOLUME_FRACTION);

#endif // ENABLE_VISUALIZE


        NuTo::CrankNicolsonEvaluate TimeIntegrationScheme(&MTStructure1D);

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


        if(measureTime)
        {
            gettimeofday(&time_end, NULL);
        }



        // %%%%%%%%%%%%%%%%%
        // Visualize Results
        // %%%%%%%%%%%%%%%%%


        if(UseVisualization)
        {
#ifdef ENABLE_VISUALIZE
            mkdir(VTKFolder.c_str(),0777);
            MTStructure1D.ExportVtkDataFileElements(VTKFolder+"/"+VTKFile,false);
#endif //ENABLE_VISUALIZE
        }

        // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        // Compare Results with paper from Johannesson and Nyman(2010)
        // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        // Just in case one uses an interpolation order higher than 1
        NNodes = NElements + 1;

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
                WPF(i)  = MTStructure1D.NodeGetNodePtr(i)->GetWaterVolumeFraction();
            }

            NuTo::FullVector<double,Eigen::Dynamic> Diff = WPF - PaperValues;
            if(std::abs(Diff.Max()) > 0.01 || std::abs(Diff.Min()) > 0.01)
            {
                throw NuTo::Exception("[Testfile: MoistureTransport1D.cpp]: Results differ to much from given values!");
            }
            else
            {
                std::cout << "Calculation finished!" << std::endl;
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
