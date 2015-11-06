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
    /** An instance of the cohort functional group definitions for this model */
    FunctionalGroupDefinitions CohortFunctionalGroupDefinitions;
    /** An instance of the stock functional group definitions for this model */
    FunctionalGroupDefinitions StockFunctionalGroupDefinitions;

    /** \brief An instance of ModelGrid to hold the grid to be used in this model */
    ModelGrid EcosystemModelGrid;

    /** \brief The number of time steps in the model run */
    unsigned NumTimeSteps;
    /** \brief The current time step */
    unsigned CurrentTimeStep;
    /** \brief The current month: 1=Jan; 2=Feb; 3=Mar etc. */
    unsigned CurrentMonth;
    /** \brief Whether to use randomisation in the model run, i.e. cohorts will be seeded with random masses and cohorts will act in a random order
     Default is true */
    bool DrawRandomly = true;
    /** \brief The threshold abundance below which cohorts will automatically become extinct */
    double ExtinctionThreshold;
    //Values to define when cohorts can be merged
    /** \brief The proportional difference in adult, juvenile and current body masses that cohorts must fall within in order to be considered for merging */
    double MergeDifference;
    /** \brief The time step units for this model */
    string GlobalModelTimeStepUnit;
    /** \brief A list of global diagnostics for this model run */
    map<string, double> GlobalDiagnosticVariables;
    /** \brief An instance of StopWatch to time individual time steps */
    StopWatch TimeStepTimer;
    StopWatch EcologyTimer;
    StopWatch OutputTimer;
    /** \brief An array of instances of the output class to deal with grid cell outputs */
    //vector<OutputCell> CellOutputs;

    /** \brief An instance of the output class to deal with global outputs */
    //OutputGlobal GlobalOutputs;
    /** \brief An instance of the output class to deal with gridded outputs */
    //OutputGrid GridOutputs;

    /** \brief A sorted list of strings from the initialisation file */
    map<string, string> InitialisationFileStrings;
    /** \brief A sorted list of strings for environmental data units */
    map<string, string> EnvironmentalDataUnits;
    /** \brief Get the human appropriation of NPP scenario to use */
    string HumanNPPExtraction;
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
    MadingleyModelInitialisation initialisation;
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
        initialisation = MadingleyModelInitialisation(
                initialisationFileName,
                OutputPath, 
                NextCohortID,
                GlobalDiagnosticVariables["NumberOfCohortsInModel"],
                GlobalDiagnosticVariables["NumberOfStocksInModel"],
                EcosystemModelGrid);

        // Assign the properties for this model run
        setUpModelRunProperties(initialisation);

        // Set up model outputs
        SetUpOutputs(initialisation);
        // Make the initial outputs
        InitialOutputs(initialisation, CurrentMonth);

        //end of initialisations
        // Initialise the cohort merger - this is just to set where the random seed comes from
        CohortMerger.SetRandom(DrawRandomly);
        // Initialise cross grid cell ecology
        disperser.setup(DrawRandomly, GlobalModelTimeStepUnit, initialisation.PlanktonDispersalThreshold);

    }
    //----------------------------------------------------------------------------------------------
    /** \brief  Run the global ecosystem model
    @param initialisation The initialization details for the current set of model simulations
     */
    void RunMadingley() {
        // Write out model run details to the console
        cout << "Running model" << endl;
        cout << "Number of time steps is: " << NumTimeSteps << endl;

        // Temporary variable
        bool varExists;
        Dispersals = 0;         
        /// Run the model
        //for (unsigned hh = 0; hh < NumTimeSteps; hh += 1) {
        for (unsigned hh = 0; hh < 2; hh += 1) {
            cout << "Running time step " << hh + 1 << "..." << endl;
            // Start the timer
            TimeStepTimer.Start();
            // Get current time step and month
            CurrentTimeStep = hh;
            CurrentMonth = Utilities.GetCurrentMonth(hh, GlobalModelTimeStepUnit);
            EcologyTimer.Start();

            RunWithinCells();

            EcologyTimer.Stop();
            cout << "Within grid ecology took: " << EcologyTimer.GetElapsedTimeSecs() << endl;

            EcologyTimer.Start();

            RunCrossGridCellEcology(Dispersals);

            EcologyTimer.Stop();
            cout << "Across grid ecology took: " << EcologyTimer.GetElapsedTimeSecs() << endl;               
            // Stop the timer
            TimeStepTimer.Stop();
            //
            OutputTimer.Start();
            // Write the global outputs for this time step
            // GlobalOutputs.TimeStepOutputs(EcosystemModelGrid, CurrentTimeStep, CurrentMonth, TimeStepTimer,CohortFunctionalGroupDefinitions,
            //                     StockFunctionalGroupDefinitions,_CellList,GlobalDiagnosticVariables, initialisation);
            OutputTimer.Stop();
            cout << "Global Outputs took: " << OutputTimer.GetElapsedTimeSecs() << endl;
            OutputTimer.Start();
            //                 {
            //                     // Write out grid outputs for this time step
            //                     GridOutputs.TimeStepOutputs(EcosystemModelGrid, CohortFunctionalGroupDefinitions, StockFunctionalGroupDefinitions, _CellList,
            //                         CurrentTimeStep, initialisation);
            //                 }
            //
            //
            OutputTimer.Stop();
            cout << "Cell/Grid Outputs took: " << OutputTimer.GetElapsedTimeSecs() << endl;
            // Write the results of dispersal to the console
            cout << "Number of cohorts that dispersed this time step: " << Dispersals << endl;
        }

        //            // Write the final global outputs
        //            GlobalOutputs.FinalOutputs();
           
    }//----------------------------------------------------------------------------------------------
    /** \brief  Run processes for cells*/
    void RunWithinCells() {
        // Instantiate a class to hold thread locked global diagnostic variables
        ThreadLockedParallelVariables singleThreadDiagnostics(0, 0, 0, NextCohortID);

        EcosystemModelGrid.ask([&](GridCell& c) {

            // Create a temporary internal copy of the grid cell stocks
            GridCellStockHandler& WorkingGridCellStocks = EcosystemModelGrid.GetGridCellStocks(c.CellEnvironment["LatIndex"][0], c.CellEnvironment["LonIndex"][0]);

            RunWithinCellStockEcology(c, WorkingGridCellStocks);
 
            RunWithinCellCohortEcology(c, singleThreadDiagnostics);
        });
        // Update the variable tracking cohort unique IDs
        NextCohortID = singleThreadDiagnostics.NextCohortIDThreadLocked;
        // Take the results from the thread local variables and apply to the global diagnostic variables
        GlobalDiagnosticVariables["NumberOfCohortsExtinct"] = singleThreadDiagnostics.Extinctions - singleThreadDiagnostics.Combinations;
        GlobalDiagnosticVariables["NumberOfCohortsProduced"] = singleThreadDiagnostics.Productions;
        GlobalDiagnosticVariables["NumberOfCohortsInModel"] = GlobalDiagnosticVariables["NumberOfCohortsInModel"] + singleThreadDiagnostics.Productions - singleThreadDiagnostics.Extinctions;
        GlobalDiagnosticVariables["NumberOfCohortsCombined"] = singleThreadDiagnostics.Combinations;
    }
    //----------------------------------------------------------------------------------------------
    /** \brief   Run ecological processes for stocks in a specified grid cell
    @param latCellIndex The latitudinal index of the cell to run stock ecology for
    @param lonCellIndex The longitudinal index of the cell to run stock ecology for
    @param workingGridCellStocks A (copy of?) the stocks in the current grid cell
     * NB - acting on a copy here? needed of parallelity? Depends if you write to the grid cell pools...
     */
    void RunWithinCellStockEcology(GridCell& gcl, GridCellStockHandler& workingGridCellStocks) {
        unsigned latCellIndex= gcl.CellEnvironment["LatIndex"][0];
        unsigned lonCellIndex= gcl.CellEnvironment["LonIndex"][0];
        // Create a local instance of the stock ecology class
        EcologyStock MadingleyEcologyStock;
        //The location of the acting stock
        vector<int> ActingStock(2);
        // Get the list of functional group indices for autotroph stocks
        vector<int> AutotrophStockFunctionalGroups = StockFunctionalGroupDefinitions.GetFunctionalGroupIndex("Heterotroph/Autotroph", "Autotroph", false);
        // Loop over autotroph functional groups
        for (int FunctionalGroup : AutotrophStockFunctionalGroups) {
            for (int ll = 0; ll < workingGridCellStocks[FunctionalGroup].size(); ll++) {
                // Get the position of the acting stock
                ActingStock[0] = FunctionalGroup;
                ActingStock[1] = ll;
                // Run stock ecology
                MadingleyEcologyStock.RunWithinCellEcology(workingGridCellStocks, ActingStock, EcosystemModelGrid.GetCellEnvironment(
                        latCellIndex, lonCellIndex), EnvironmentalDataUnits, HumanNPPExtraction, StockFunctionalGroupDefinitions,
                        CurrentTimeStep, GlobalModelTimeStepUnit,CurrentMonth,
                        initialisation.InitialisationFileStrings["OutputDetail"]);
                //workingGridCellStocks[ActingStock].TotalBiomass *= 0.75;//MB strange line - commented out in original?
            }
        }

    }
    //----------------------------------------------------------------------------------------------
    /** \brief   Run ecological processes for cohorts in a specified grid cell
    @param gcl Reference to the current grid cell
    @param partial Trakc some global variables pertaining to cohort numbers etc.

     * NB - need to take care here when cohort updates get applied

    */
    void RunWithinCellCohortEcology(GridCell& gcl, ThreadLockedParallelVariables& partial) {
        // Local instances of classes
        EcologyCohort MadingleyEcologyCohort;

        // Initialize ecology for stocks and cohorts
        MadingleyEcologyCohort.InitializeEcology(gcl.CellEnvironment["Cell Area"][0],GlobalModelTimeStepUnit, DrawRandomly);
        Activity CohortActivity;
        
        // Diagnostic biological variables don't need to be reset every cohort, but rather every grid cell
        int EcosystemModelParallelTempval1=0,EcosystemModelParallelTempval2 = 0;

        // Initialise eating formulations - has to be redone ecery step?
        MadingleyEcologyCohort.EatingFormulations["Basic eating"]->InitializeEcologicalProcess(gcl.GridCellCohorts, gcl.GridCellStocks,
                CohortFunctionalGroupDefinitions, StockFunctionalGroupDefinitions, "revised predation");

        MadingleyEcologyCohort.EatingFormulations["Basic eating"]->InitializeEcologicalProcess(gcl.GridCellCohorts, gcl.GridCellStocks
                , CohortFunctionalGroupDefinitions, StockFunctionalGroupDefinitions, "revised herbivory");

        // Loop over randomly ordered gridCellCohorts to implement biological functions
        //if (DrawRandomly) {
        // Randomly order the cohort indices
        //    RandomCohortOrder = Utilities.RandomlyOrderedCohorts(TotalCohortNumber);
        //} else {
        //    RandomCohortOrder = Utilities.NonRandomlyOrderedCohorts(TotalCohortNumber, CurrentTimeStep);
        //}
        //MB need to put back random ordering
        gcl.ask([&](Cohort& c){

            // Perform all biological functions except dispersal (which is cross grid cell)
            if (gcl.GridCellCohorts[c.FunctionalGroupIndex].size() != 0 && c.CohortAbundance > ExtinctionThreshold) {
                // Calculate number of cohorts in this functional group in this grid cell before running ecology
                EcosystemModelParallelTempval1 = gcl.GridCellCohorts[c.FunctionalGroupIndex].size();

                CohortActivity.AssignProportionTimeActive(c, gcl.CellEnvironment, CohortFunctionalGroupDefinitions, CurrentTimeStep, CurrentMonth);

                // Run ecology
                MadingleyEcologyCohort.RunWithinCellEcology(gcl,c,
                        CohortFunctionalGroupDefinitions, StockFunctionalGroupDefinitions, CurrentTimeStep,
                         partial,  CurrentMonth, initialisation);

                // Update the properties of the acting cohort
                MadingleyEcologyCohort.UpdateEcology(gcl, c,
                         CohortFunctionalGroupDefinitions, StockFunctionalGroupDefinitions, CurrentTimeStep);

                // Add newly produced cohorts to the tracking variable
                EcosystemModelParallelTempval2 += gcl.GridCellCohorts[c.FunctionalGroupIndex].size() - EcosystemModelParallelTempval1;


                // Check that the mass of individuals in this cohort is still >= 0 after running ecology
                assert(c.IndividualBodyMass >= 0.0 && "Biomass < 0 for this cohort");
            }

            // Check that the mass of individuals in this cohort is still >= 0 after running ecology
            if (gcl.GridCellCohorts[c.FunctionalGroupIndex].size() > 0)assert(c.IndividualBodyMass >= 0.0 && "Biomass < 0 for this cohort");
        });


        // Update diagnostics of productions
        partial.Productions += EcosystemModelParallelTempval2;

        RunExtinction(gcl, partial);

        // Merge cohorts, if necessary
        if (gcl.GridCellCohorts.GetNumberOfCohorts() > initialisation.MaxNumberOfCohorts) {
            partial.Combinations = CohortMerger.MergeToReachThresholdFast(gcl.GridCellCohorts, gcl.GridCellCohorts.GetNumberOfCohorts(), initialisation.MaxNumberOfCohorts);

            //Run extinction a second time to remove those cohorts that have been set to zero abundance when merging
            RunExtinction(gcl, partial);
        } else
        partial.Combinations = 0;

        // Write out the updated cohort numbers after all ecological processes have occurred
        //MB this is changed!
        //EcosystemModelGrid.SetGridCellCohorts(gcl.GridCellCohorts, latCellIndex, lonCellIndex);
    }

    //----------------------------------------------------------------------------------------------
    /** \brief Carries out extinction on cohorts that have an abundance below a defined extinction threshold */
    void RunExtinction(GridCell& gcl, ThreadLockedParallelVariables& partial) {
        bool VarExists;

        // Loop over cohorts and remove any whose abundance is below the extinction threshold
        for (int FG = 0; FG < CohortFunctionalGroupDefinitions.GetNumberOfFunctionalGroups(); FG++) {
            // Create a list to hold the cohorts to remove
            vector<Cohort>CohortsToRemove;
            // Loop through each cohort in the functional group
            for(auto& c:gcl.GridCellCohorts[FG] ){
                // If cohort abundance is less than the extinction threshold then add to the list for extinction
                if (c.CohortAbundance < ExtinctionThreshold || c.IndividualBodyMass <= 1.e-300) {
                    CohortsToRemove.push_back(c);

                    partial.Extinctions += 1;

                }
            }
            //cohort order needs to be reversed to get removals right
            reverse(CohortsToRemove.begin(),CohortsToRemove.end());
            // Code to add the biomass to the biomass pool and dispose of the cohort
            for (auto& c :CohortsToRemove) {
                // Add biomass of the extinct cohort to the organic matter pool
                unsigned latCellIndex=gcl.CellEnvironment["LatIndex"][0],lonCellIndex=gcl.CellEnvironment["LonIndex"][0];
                EcosystemModelGrid.SetEnviroLayer("Organic Pool", 0, EcosystemModelGrid.GetEnviroLayer("Organic Pool", 0, latCellIndex, lonCellIndex, VarExists) +
                        (c.IndividualBodyMass + c.IndividualReproductivePotentialMass) * c.CohortAbundance, latCellIndex, lonCellIndex);
                assert(EcosystemModelGrid.GetEnviroLayer("Organic Pool", 0, latCellIndex, lonCellIndex, VarExists) >= 0 && "Organic pool < 0");
 
                // Remove the extinct cohort from the list of cohorts
                EcosystemModelGrid.DeleteGridCellIndividualCohort(c);

            }

        }
    }
    //----------------------------------------------------------------------------------------------
    /** \brief Run ecological processes that operate across grid cells */
    void RunCrossGridCellEcology(unsigned& dispersals) {
        // Loop through each grid cell, and run dispersal for each.
        // Note that currently dispersal is not parallelised, although it could be (though care would need to be taken to ensure that necessary variables are thread-locked
        
        EcosystemModelGrid.ask([&](GridCell& c) {

            // We have looped through individual cells and calculated ecological processes for each. Now do this for cross grid cell dispersal

            disperser.RunCrossGridCellEcologicalProcess(c, EcosystemModelGrid,  CohortFunctionalGroupDefinitions,  CurrentMonth);

        });
        // Apply the changes from dispersal
        disperser.UpdateCrossGridCellEcology(EcosystemModelGrid, dispersals);

    }
    //----------------------------------------------------------------------------------------------
    /** \brief  Assigns the properties of the current model run
    @param initialisation An instance of the model initialisation class 
    @param scenarioParameters The parameters for the scenarios to run
    @param scenarioIndex The index of the scenario that this model is to run
    @param outputFilesSuffix The suffix to be applied to all outputs from this model run
     */
    void setUpModelRunProperties(MadingleyModelInitialisation& initialisation) {
        // Assign the properties of this model run from the same properties in the specified model initialisation
        GlobalModelTimeStepUnit = initialisation.GlobalModelTimeStepUnit;
        NumTimeSteps = initialisation.NumTimeSteps;

        DrawRandomly = initialisation.DrawRandomly;
        ExtinctionThreshold = initialisation.ExtinctionThreshold;
        MergeDifference = initialisation.MergeDifference;
        InitialisationFileStrings = initialisation.InitialisationFileStrings;
        CohortFunctionalGroupDefinitions = initialisation.CohortFunctionalGroupDefinitions;
        StockFunctionalGroupDefinitions = initialisation.StockFunctionalGroupDefinitions;

        HumanNPPExtraction = initialisation.InitialisationFileStrings["HumanNPPExtraction"];

        EnvironmentalDataUnits = initialisation.Units;

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
    void SetUpOutputs(MadingleyModelInitialisation initialisation) {
        //            // Initialise the global outputs
        //            GlobalOutputs = new OutputGlobal(InitialisationFileStrings["OutputDetail"], initialisation);
        //

        //                GridOutputs = new OutputGrid(InitialisationFileStrings["OutputDetail"], initialisation);
        //            
        //
        //            
    }
    //----------------------------------------------------------------------------------------------
    /** \brief   Generates the initial outputs for this model run
    @param outputFilesSuffix The suffix to be applied to all outputs from this model run
     */
    void InitialOutputs(MadingleyModelInitialisation initialisation, unsigned month) {
        //            // Set up global outputs for all model runs
        //            GlobalOutputs.SetupOutputs(NumTimeSteps, EcosystemModelGrid, OutputFilesSuffix);
        //
        //            // Create initial global outputs
        //            GlobalOutputs.InitialOutputs(EcosystemModelGrid,CohortFunctionalGroupDefinitions,StockFunctionalGroupDefinitions,_CellList,
        //                GlobalDiagnosticVariables, initialisation);
        //
        //            // Temporary
        //            Boolean varExists;
        //

        //                // Set up grid outputs
        //                GridOutputs.SetupOutputs(EcosystemModelGrid, OutputFilesSuffix, NumTimeSteps, 
        //                    CohortFunctionalGroupDefinitions, StockFunctionalGroupDefinitions);
        //
        //                // Create initial grid outputs
        //                GridOutputs.InitialOutputs(EcosystemModelGrid, CohortFunctionalGroupDefinitions, StockFunctionalGroupDefinitions, _CellList, initialisation);
        //            
        //
    }
    

};
#endif