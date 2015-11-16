#ifndef MADINGLEYMODEL_H
#define MADINGLEYMODEL_H
#include <vector>
#include <map>
#include <MadingleyModelInitialisation.h>
#include <FunctionalGroupDefinitions.h>
#include <Stopwatch.h>
#include <CohortMerge.h>
#include <ModelGrid.h>
#include <GridCell.h>
#include <Dispersal.h>
#include <EcologyStock.h>
#include <EcologyCohort.h>
#include <Activity.h>
#include <ThreadLocked.h>
#include <fstream>
/// @todo check private versus public variables
/** \file MadingleyModel.h
 * \brief The main model header file
 * */

/** \brief The ecosystem model */
class MadingleyModel {
public:
    //----------------------------------------------------------------------------------------------
    //Variables
    //----------------------------------------------------------------------------------------------       

    /** \brief An instance of ModelGrid to hold the grid to be used in this model */
    ModelGrid EcosystemModelGrid;

    /** \brief The current time step */
    unsigned CurrentTimeStep;
    /** \brief The current month: 1=Jan; 2=Feb; 3=Mar etc. */
    unsigned CurrentMonth;

    /** \brief A list of global diagnostics for this model run */
    map<string, double> GlobalDiagnosticVariables;
    /** \brief An instance of StopWatch to time individual time steps */
    StopWatch TimeStepTimer;
    StopWatch EcologyTimer;
    StopWatch DispersalTimer;
    StopWatch OutputTimer;
    StopWatch MergeTimer;
    /** \brief An array of instances of the output class to deal with grid cell outputs */
    //vector<OutputCell> CellOutputs;

    /** \brief An instance of the output class to deal with global outputs */
    //OutputGlobal GlobalOutputs;
    /** \brief An instance of the output class to deal with gridded outputs */
    //OutputGrid GridOutputs;

    /** A variable to increment for the purposes of giving each cohort a unique ID */
    long long NextCohortID; 
    /** \brief Variable to track the number of cohorts that have dispersed. Doesn't need to be thread-local because all threads have converged prior to running cross-grid-cell processes */
    unsigned Dispersals;
    /** \brief Instance of the class to perform general functions */
    UtilityFunctions Utilities;
    /** \brief An instance of the merging class */
    CohortMerge CohortMerger;
    /** \brief An instance of the simple random number generator class */
    std::default_random_engine RandomNumberGenerator;
    //
    MadingleyModelInitialisation params;
    Dispersal disperser;
    //----------------------------------------------------------------------------------------------
    //Methods
    //----------------------------------------------------------------------------------------------
    
    //----------------------------------------------------------------------------------------------
    /** \brief   Initializes the ecosystem model
    @param initialisationFileName The name of the file with model parameters
    @param OutputPath Where the output will be stored
     */
    MadingleyModel(string initialisationFileName, string OutputPath) {
        // Set up list of global diagnostics
        SetUpGlobalDiagnosticsList();
        // Initialise the cohort ID to zero
        NextCohortID = 0;
        params = MadingleyModelInitialisation(
                initialisationFileName,
                OutputPath, 
                NextCohortID,
                GlobalDiagnosticVariables["NumberOfCohortsInModel"],
                GlobalDiagnosticVariables["NumberOfStocksInModel"],
                EcosystemModelGrid);


        // Set up model outputs
        SetUpOutputs();
        // Make the initial outputs
        InitialOutputs( CurrentMonth);

        //end of initialisations
        // Initialise the cohort merger - this is just to set where the random seed comes from
        CohortMerger.SetRandom(params.DrawRandomly);
        // Initialise cross grid cell ecology
        disperser.setup(params);

    }
    //----------------------------------------------------------------------------------------------
    /** \brief  Run the global ecosystem model     */
    void RunMadingley() {
        // Write out model run details to the console
        cout << "Running model" << endl;
        cout << "Number of time steps is: " << params.NumTimeSteps << endl;
        ofstream urg("aaggh");
        urg<<"step,dispersals,extinctions,productions,combinations,totalcohorts,incelltime,dispersaltime"<<endl;
        // Temporary variable
        bool varExists;
        Dispersals = 0;         
        /// Run the model
        for (unsigned hh = 0; hh < params.NumTimeSteps; hh += 1) {
        //for (unsigned hh = 0; hh < 2; hh += 1) {
            cout << "Running time step " << hh + 1 << "..." << endl;
            // Start the timer
            TimeStepTimer.Start();
            // Get current time step and month
            CurrentTimeStep = hh;
            CurrentMonth = Utilities.GetCurrentMonth(hh, params.GlobalModelTimeStepUnit);
            EcologyTimer.Start();

            RunWithinCells();

            EcologyTimer.Stop();
            cout << "Within grid ecology took: " << EcologyTimer.GetElapsedTimeSecs() << endl;

            DispersalTimer.Start();

            RunCrossGridCellEcology(Dispersals);

            DispersalTimer.Stop();
            cout << "Across grid ecology took: " << DispersalTimer.GetElapsedTimeSecs() << endl;               
            // Stop the timer
            TimeStepTimer.Stop();
            //
            OutputTimer.Start();
            // Write the global outputs for this time step
            // GlobalOutputs.TimeStepOutputs(EcosystemModelGrid, CurrentTimeStep, CurrentMonth, TimeStepTimer,CohortFunctionalGroupDefinitions,
            //                     StockFunctionalGroupDefinitions,_CellList,GlobalDiagnosticVariables, params);
            OutputTimer.Stop();
            cout << "Global Outputs took: " << OutputTimer.GetElapsedTimeSecs() << endl;
            OutputTimer.Start();
            //                 {
            //                     // Write out grid outputs for this time step
            //                     GridOutputs.TimeStepOutputs(EcosystemModelGrid, CohortFunctionalGroupDefinitions, StockFunctionalGroupDefinitions, _CellList,
            //                         CurrentTimeStep, params);
            //                 }
            //
            //
            OutputTimer.Stop();
            cout << "Cell/Grid Outputs took: " << OutputTimer.GetElapsedTimeSecs() << endl;
            // Write the results of dispersal to the console

        cout<<"Total Cohorts remaining "<<GlobalDiagnosticVariables["NumberOfCohortsInModel"]<<endl ;
        urg<<hh<<","<<Dispersals<<","<<GlobalDiagnosticVariables["NumberOfCohortsExtinct"]<<","<<GlobalDiagnosticVariables["NumberOfCohortsProduced"]<<","<<GlobalDiagnosticVariables["NumberOfCohortsCombined"]<<","<<GlobalDiagnosticVariables["NumberOfCohortsInModel"]<<","<<EcologyTimer.GetElapsedTimeSecs()<<","<<DispersalTimer.GetElapsedTimeSecs()<<endl;

        }

        //            // Write the final global outputs
        //            GlobalOutputs.FinalOutputs();
           
    }//----------------------------------------------------------------------------------------------
    /** \brief  Run processes for cells*/
    void RunWithinCells() {
        // Instantiate a class to hold thread locked global diagnostic variables
        ThreadLockedParallelVariables singleThreadDiagnostics(0, 0, 0, NextCohortID);

        EcosystemModelGrid.ask([&](GridCell& gcl) {

            RunWithinCellStockEcology(gcl);
 
            RunWithinCellCohortEcology(gcl, singleThreadDiagnostics);
        });
        // Update the variable tracking cohort unique IDs
        NextCohortID = singleThreadDiagnostics.NextCohortIDThreadLocked;
        // Take the results from the thread local variables and apply to the global diagnostic variables
        GlobalDiagnosticVariables["NumberOfCohortsExtinct"] = singleThreadDiagnostics.Extinctions - singleThreadDiagnostics.Combinations;
        GlobalDiagnosticVariables["NumberOfCohortsProduced"] = singleThreadDiagnostics.Productions;
        GlobalDiagnosticVariables["NumberOfCohortsInModel"] = GlobalDiagnosticVariables["NumberOfCohortsInModel"] + singleThreadDiagnostics.Productions - singleThreadDiagnostics.Extinctions;
        GlobalDiagnosticVariables["NumberOfCohortsCombined"] = singleThreadDiagnostics.Combinations;
        cout << "Number of cohorts that dispersed this time step: " << Dispersals << endl;
        cout<<"Extinctions this step "<<singleThreadDiagnostics.Extinctions - singleThreadDiagnostics.Combinations<<endl ;
        cout<<"Productions "<<singleThreadDiagnostics.Productions<<endl ;
        cout<<"Combinations"<<singleThreadDiagnostics.Combinations<<endl ;
    }
    //----------------------------------------------------------------------------------------------
    /** \brief   Run ecological processes for stocks in a specified grid cell
    @param gcl The current cell 
     */
    void RunWithinCellStockEcology(GridCell& gcl) {

        // Create a local instance of the stock ecology class
        EcologyStock MadingleyEcologyStock;
        // Get the list of functional group indices for autotroph stocks
        vector<int> AutotrophStockFunctionalGroups = params.StockFunctionalGroupDefinitions.GetFunctionalGroupIndex("Heterotroph/Autotroph", "Autotroph", false);
        // Loop over autotroph functional groups
        for (unsigned FunctionalGroup : AutotrophStockFunctionalGroups) {
            for (auto& ActingStock: gcl.GridCellStocks[FunctionalGroup]) {

                // Run stock ecology
                MadingleyEcologyStock.RunWithinCellEcology(gcl,ActingStock,
                        CurrentTimeStep, CurrentMonth,params);
                //workingGridCellStocks[ActingStock].TotalBiomass *= 0.75;//MB strange line - commented out in original?
            }
        }

    }
    //----------------------------------------------------------------------------------------------
    /** \brief   Run ecological processes for cohorts in a specified grid cell
    @param gcl Reference to the current grid cell
    @param partial Track some global variables pertaining to cohort numbers etc.

     * NB - need to take care here when cohort updates get applied

    */
    void RunWithinCellCohortEcology(GridCell& gcl, ThreadLockedParallelVariables& partial) {
        // Local instances of classes
        // Initialize ecology for stocks and cohorts - needed fresh every timestep?

        EcologyCohort MadingleyEcologyCohort(gcl,params);

        Activity CohortActivity;
        
        // Loop over randomly ordered gridCellCohorts to implement biological functions
        //if (DrawRandomly) {
        // Randomly order the cohort indices
        //    RandomCohortOrder = Utilities.RandomlyOrderedCohorts(TotalCohortNumber);
        //} else {
        //    RandomCohortOrder = Utilities.NonRandomlyOrderedCohorts(TotalCohortNumber, CurrentTimeStep);
        //}
        //MB need to put back random ordering - NB apply ecology may set bodymass of a cohort to zero!
        gcl.ask([&](Cohort& c){
            // Perform all biological functions except dispersal (which is cross grid cell)
            if (gcl.GridCellCohorts[c.FunctionalGroupIndex].size() != 0 && c.CohortAbundance > params.ExtinctionThreshold) {

                CohortActivity.AssignProportionTimeActive(gcl,c, CurrentTimeStep, CurrentMonth, params);

                // Run ecology
                MadingleyEcologyCohort.RunWithinCellEcology(gcl,c,CurrentTimeStep, partial,  CurrentMonth, params);

                // Update the properties of the acting cohort
                MadingleyEcologyCohort.UpdateEcology(gcl, c, CurrentTimeStep);
                Cohort::zeroDeltas();
                
                // Check that the mass of individuals in this cohort is still >= 0 after running ecology
                assert(c.IndividualBodyMass >= 0.0 && "Biomass < 0 for this cohort");
            }

            // Check that the mass of individuals in this cohort is still >= 0 after running ecology
            if (gcl.GridCellCohorts[c.FunctionalGroupIndex].size() > 0)assert(c.IndividualBodyMass >= 0.0 && "Biomass < 0 for this cohort");
        });
        

        for (auto& c: Cohort::newCohorts){
            EcosystemModelGrid.AddNewCohortToGridCell(c);
            if (c.destination != &gcl)cout<<"whut? wrong cell?"<<endl;
        }
        partial.Productions+=Cohort::newCohorts.size();
        Cohort::newCohorts.clear();

        RunExtinction(gcl, partial);

        // Merge cohorts, if necessary
        if (gcl.GetNumberOfCohorts() > params.MaxNumberOfCohorts) {
            partial.Combinations += CohortMerger.MergeToReachThresholdFast(gcl, params);
         

            //Run extinction a second time to remove those cohorts that have been set to zero abundance when merging
            RunExtinction(gcl, partial);

        }
    }

    //----------------------------------------------------------------------------------------------
    /** \brief Carries out extinction on cohorts that have an abundance below a defined extinction threshold */
    void RunExtinction(GridCell& gcl, ThreadLockedParallelVariables& partial) {
        bool VarExists;
        // Loop over cohorts and remove any whose abundance is below the extinction threshold

        vector<Cohort>CohortsToRemove;
        gcl.ask([&](Cohort & c) {
            if (c.CohortAbundance < params.ExtinctionThreshold || c.IndividualBodyMass <= 1.e-300) {
                CohortsToRemove.push_back(c);
                partial.Extinctions += 1;}
            });

            // Code to add the biomass to the biomass pool and dispose of the cohort
            for (auto& c :CohortsToRemove) {
                // Add biomass of the extinct cohort to the organic matter pool

                double deadMatter=(c.IndividualBodyMass + c.IndividualReproductivePotentialMass) * c.CohortAbundance;
                EcosystemModelGrid.AddToEnviroLayer("Organic Pool", 0, deadMatter , c);
                assert(EcosystemModelGrid.GetEnviroLayer("Organic Pool", 0, c, VarExists) >= 0 && "Organic pool < 0");
 
                // Remove the extinct cohort from the list of cohorts
                EcosystemModelGrid.DeleteGridCellIndividualCohort(c);

            }

    }
    //----------------------------------------------------------------------------------------------
    /** \brief Run ecological processes that operate across grid cells */
    void RunCrossGridCellEcology(unsigned& dispersals) {
        // Loop through each grid cell, and run dispersal for each.

        
        EcosystemModelGrid.ask([&](GridCell& c) {

            disperser.RunCrossGridCellEcologicalProcess(c, EcosystemModelGrid,  params,  CurrentMonth);

        });
        
        // Apply the changes from dispersal
        disperser.UpdateCrossGridCellEcology(EcosystemModelGrid, dispersals);
    }
    //----------------------------------------------------------------------------------------------
    /** \brief   Sets up the list of global diagnostic variables
     */
    void SetUpGlobalDiagnosticsList() {
        // Add global diagnostic variables
        GlobalDiagnosticVariables["NumberOfCohortsExtinct"] = 0.0;
        GlobalDiagnosticVariables["NumberOfCohortsProduced"] = 0.0;
        GlobalDiagnosticVariables["NumberOfCohortsCombined"] = 0.0;
        GlobalDiagnosticVariables["NumberOfCohortsInModel"] = 0.0;
        GlobalDiagnosticVariables["NumberOfStocksInModel"] = 0.0;
    }
    //----------------------------------------------------------------------------------------------
    /** \brief   Sets up the model outputs
    @param initialisation An instance of the model initialisation class
     */
    void SetUpOutputs() {
        //            // Initialise the global outputs
        //            GlobalOutputs = new OutputGlobal(InitialisationFileStrings["OutputDetail"], params);
        //

        //                GridOutputs = new OutputGrid(InitialisationFileStrings["OutputDetail"], params);
        //            
        //
        //            
    }
    //----------------------------------------------------------------------------------------------
    /** \brief   Generates the initial outputs for this model run
    @param outputFilesSuffix The suffix to be applied to all outputs from this model run
     */
    void InitialOutputs( unsigned month) {
        //            // Set up global outputs for all model runs
        //            GlobalOutputs.SetupOutputs(NumTimeSteps, EcosystemModelGrid, OutputFilesSuffix);
        //
        //            // Create initial global outputs
        //            GlobalOutputs.InitialOutputs(EcosystemModelGrid,CohortFunctionalGroupDefinitions,StockFunctionalGroupDefinitions,_CellList,
        //                GlobalDiagnosticVariables, params);
        //
        //            // Temporary
        //            Boolean varExists;
        //

        //                // Set up grid outputs
        //                GridOutputs.SetupOutputs(EcosystemModelGrid, OutputFilesSuffix, NumTimeSteps, 
        //                    CohortFunctionalGroupDefinitions, StockFunctionalGroupDefinitions);
        //
        //                // Create initial grid outputs
        //                GridOutputs.InitialOutputs(EcosystemModelGrid, CohortFunctionalGroupDefinitions, StockFunctionalGroupDefinitions, _CellList, params);
        //            
        //
    }
    

};
#endif