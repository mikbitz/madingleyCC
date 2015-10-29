#ifndef MADINGLEYMODEL_H
#define MADINGLEYMODEL_H
#include <vector>
#include <map>
#include <MadingleyModelInitialisation.h>
#include <FunctionalGroupDefinitions.h>
#include <Stopwatch.h>
#include <CohortMerge.h>
#include <ModelGrid.h>
#include <Dispersal.h>
#include <EcologyStock.h>
#include <EcologyCohort.h>
#include <EnviroData.h>
#include <ProcessTracker.h>
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
    /** \brief A list of environmental data layers */
    map<string, EnviroData> EnviroStack;
    /** \brief An instance of ModelGrid to hold the grid to be used in this model */
    ModelGrid EcosystemModelGrid;
    /** The lowest latitude for the model grid */
    float BottomLatitude;
    /** The upper latitude for the model grid */
    float TopLatitude;
    /**  The left-most longitude for the model grid */
    float LeftmostLongitude;
    /**  The right-most longitude for the model grid */
    float RightmostLongitude;
    /** The size of the grid cells in degrees */
    float CellSize;
    /** \brief The rarefaction to be applied to live cells in the model grid */
    int CellRarefaction;
    /** \brief The number of time steps in the model run */
    unsigned NumTimeSteps;
    /** \brief The current time step */
    unsigned CurrentTimeStep;
    /** \brief The current month: 1=Jan; 2=Feb; 3=Mar etc. */
    unsigned CurrentMonth;
    /** \brief Whether to use randomisation in the model run, i.e. cohorts will be seeeded with random masses and cohorts will act in a random order
     Default is true */
    bool DrawRandomly = true;
    /** \brief The threshold abundance below which cohorts will automatically become extinct */
    double ExtinctionThreshold;
    //Values to define when cohorts can be merged
    /** \brief The proportional difference in adult, juvenile and current body masses that cohorts must fall within in order to be considered for merging */
    double MergeDifference;
    /** \brief The time step units for this model */
    string GlobalModelTimeStepUnit;
    /** \brief Pairs of longitude and latitude indices for all active cells in the model grid */
    vector< vector<unsigned> > CellList;
    /** \brief A list of global diagnostics for this model run */
    map<string, double> GlobalDiagnosticVariables;
    /** \brief Whether the model will run in parallel (default  is false) */
    bool RunGridCellsInParallel = false;
    /** \brief An instance of StopWatch to time individual time steps */
    StopWatch TimeStepTimer;
    StopWatch EcologyTimer;
    StopWatch OutputTimer;
    /** \brief An array of instances of the output class to deal with grid cell outputs */
    //vector<OutputCell> CellOutputs;
    /** \brief  An array of indices of process trackers for each grid cell */
    vector<ProcessTracker> ProcessTrackers;
    /** \brief  An instance of a global process tracker to track across the model grid */
    GlobalProcessTracker TrackGlobalProcesses;
    /** \brief An instance of the output class to deal with global outputs */
    //OutputGlobal GlobalOutputs;
    /** \brief An instance of the output class to deal with gridded outputs */
    //OutputGrid GridOutputs;
    /** \brief The suffix to be applied to files output by this model instance */
    string OutputFilesSuffix;
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
    //
    MadingleyModelInitialisation initialisation;
    Dispersal disperser;
    //----------------------------------------------------------------------------------------------
    //Methods
    //----------------------------------------------------------------------------------------------
    
    //----------------------------------------------------------------------------------------------
    /** \brief   Initializes the ecosystem model
    @param initialisation An instance of the model initialisation class 
    @param scenarioParameters The parameters for the scenarios to run
    @param scenarioIndex The index of the scenario that this model is to run
    @param outputFilesSuffix The suffix to be applied to all outputs from this model run
    @param globalModelTimeStepUnit The time step unit used in the model
     */
    MadingleyModel(string initialisationFileName, string OutputPath) {
        initialisation = MadingleyModelInitialisation(initialisationFileName, OutputPath);
        // Assign the properties for this model run
        AssignModelRunProperties(initialisation);
        // Set up list of global diagnostics
        SetUpGlobalDiagnosticsList();
        // Set up the model grid
        SetUpModelGrid(initialisation);
        // Set up model outputs
        SetUpOutputs(initialisation);
        // Make the initial outputs
        InitialOutputs(initialisation, CurrentMonth);
        //Instance the array of process trackers
        //ProcessTrackers.resize[CellList.size()];
        // Temporary variables
        bool varExists;
        //
        //            // Set up process trackers for each grid cell
        //            for (int i = 0; i < _CellList.Count; i++)
        //            {
        //                ProcessTrackers[i] = new ProcessTracker(NumTimeSteps,
        //                EcosystemModelGrid.Lats, EcosystemModelGrid.Lons,
        //                _CellList,
        //                initialisation.ProcessTrackingOutputs,
        //                initialisation.TrackProcesses,
        //                CohortFunctionalGroupDefinitions,
        //                EcosystemModelGrid.GlobalMissingValue,
        //                outputFilesSuffix,
        //                initialisation.OutputPath, initialisation.ModelMassBins,
        //                SpecificLocations, i, initialisation, 
        //                EcosystemModelGrid.GetEnviroLayer("Realm", 0, _CellList[i][0], _CellList[i][1], out varExists) == 2.0,
        //                EcosystemModelGrid.LatCellSize,
        //                EcosystemModelGrid.LonCellSize);
        //            }
        //            
        //            //Set up a global process tracker
        //            if (SpecificLocations) initialisation.TrackGlobalProcesses = false;
        //
        //            TrackGlobalProcesses = new GlobalProcessTracker(NumTimeSteps,
        //                EcosystemModelGrid.Lats, EcosystemModelGrid.Lons,
        //                _CellList,
        //                initialisation.ProcessTrackingOutputs,
        //                initialisation.TrackGlobalProcesses,
        //                CohortFunctionalGroupDefinitions,
        //                EcosystemModelGrid.GlobalMissingValue,
        //                outputFilesSuffix,
        //                initialisation.OutputPath, initialisation.ModelMassBins,
        //                SpecificLocations, initialisation,
        //                EcosystemModelGrid.LatCellSize,
        //                EcosystemModelGrid.LonCellSize);
        //
        //Record the initial cohorts in the process trackers
        //RecordInitialCohorts();
        // Set the global model time step unit
        GlobalModelTimeStepUnit = initialisation.GlobalModelTimeStepUnit;
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
        for (unsigned hh = 0; hh < 1; hh += 1) {
            cout << "Running time step " << hh + 1 << "..." << endl;
            // Start the timer
            TimeStepTimer.Start();
            // Get current time step and month
            CurrentTimeStep = hh;
            CurrentMonth = Utilities.GetCurrentMonth(hh, GlobalModelTimeStepUnit);
            EcologyTimer.Start();
            // Run cells in sequence
            RunCellsSequentially();
            EcologyTimer.Stop();
            cout << "Within grid ecology took: " << EcologyTimer.GetElapsedTimeSecs() << endl;
            //if(TrackGlobalProcesses.TrackProcesses) TrackGlobalProcesses.StoreNPPGrid(hh);
            EcologyTimer.Start();
            // Run cross grid cell ecology
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
        //
        //             if (TrackGlobalProcesses.TrackProcesses) TrackGlobalProcesses.CloseNPPFile();
        //
        //            // Loop over cells and close process trackers
        //             for (int i = 0; i < _CellList.Count; i++)
        //             {
        //                 if (ProcessTrackers[i].TrackProcesses) ProcessTrackers[i].CloseStreams(SpecificLocations);
        //             }
        //
        //            // Write the final global outputs
        //            GlobalOutputs.FinalOutputs();
           
    }
    //----------------------------------------------------------------------------------------------
    /** \brief  A method to run the main ecosystem model loop in parallel (latitudinal strips)
    @param latCellIndex The latitude index of the cell to run
    @param lonCellIndex The longitude index of the cell to run
    @param partial A threadlockedparallelvariable that is used to pass global diagnostic information back with locking or race conditions
    @param dispersalOnly Whether to run dispersal only (i.e. to turn all other ecological processes off
    @remarks Note that variables and instances of classes that are written to within this method MUST be local within this method to prevent 
     race issues and multiple threads attempting to write to the same variable when running the program in parallel
     */
    void RunCell(int cellIndex, ThreadLockedParallelVariables& partial,  MadingleyModelInitialisation& initialisation) {
        // Create a temporary internal copy of the grid cell cohorts
        GridCellCohortHandler WorkingGridCellCohorts = EcosystemModelGrid.GetGridCellCohorts(CellList[cellIndex][0], CellList[cellIndex][1]);
        // Create a temporary internal copy of the grid cell stocks
        GridCellStockHandler WorkingGridCellStocks = EcosystemModelGrid.GetGridCellStocks(CellList[cellIndex][0], CellList[cellIndex][1]);
        // Run stock ecology
        //RunWithinCellStockEcology(CellList[cellIndex][0], CellList[cellIndex][1], WorkingGridCellStocks, cellIndex);
        // Run cohort ecology
        RunWithinCellCohortEcology(CellList[cellIndex][0], CellList[cellIndex][1], partial, WorkingGridCellCohorts, WorkingGridCellStocks, InitialisationFileStrings["OutputDetail"], cellIndex, initialisation);
    }
    //----------------------------------------------------------------------------------------------
    /** \brief  Assigns the properties of the current model run
    @param initialisation An instance of the model initialisation class 
    @param scenarioParameters The parameters for the scenarios to run
    @param scenarioIndex The index of the scenario that this model is to run
    @param outputFilesSuffix The suffix to be applied to all outputs from this model run
     */
    void AssignModelRunProperties(MadingleyModelInitialisation& initialisation) {
        // Assign the properties of this model run from the same properties in the specified model initialisation
        GlobalModelTimeStepUnit = initialisation.GlobalModelTimeStepUnit;
        NumTimeSteps = initialisation.NumTimeSteps;
        CellSize = initialisation.CellSize;
        BottomLatitude = initialisation.BottomLatitude;
        TopLatitude = initialisation.TopLatitude;
        LeftmostLongitude = initialisation.LeftmostLongitude;
        RightmostLongitude = initialisation.RightmostLongitude;
        CellRarefaction = initialisation.CellRarefaction;
        RunGridCellsInParallel = initialisation.RunCellsInParallel;
        DrawRandomly = initialisation.DrawRandomly;
        ExtinctionThreshold = initialisation.ExtinctionThreshold;
        MergeDifference = initialisation.MergeDifference;
        InitialisationFileStrings = initialisation.InitialisationFileStrings;
        CohortFunctionalGroupDefinitions = initialisation.CohortFunctionalGroupDefinitions;
        StockFunctionalGroupDefinitions = initialisation.StockFunctionalGroupDefinitions;
        EnviroStack = initialisation.EnviroStack;
        HumanNPPExtraction = initialisation.InitialisationFileStrings["HumanNPPExtraction"];
        OutputFilesSuffix = "";
        EnvironmentalDataUnits = initialisation.Units;

        // If the model run is for a whole grid, then create envirodata classes with the environmental data
            // TEMP COMMENT OUT OF TEMPERATURE
            //MB originally these were done with "FetchClimate" but this needs to be replaced
            if (EnviroStack.count("Temperature")   == 0) EnviroStack["Temperature"]   = EnviroData("temperature", "month", BottomLatitude, LeftmostLongitude, TopLatitude, RightmostLongitude, CellSize); //, EnvironmentalDataSource.ANY));
            if (EnviroStack.count("LandDTR")       == 0) EnviroStack["LandDTR"]       = EnviroData("land_dtr", "month", BottomLatitude, LeftmostLongitude, TopLatitude, RightmostLongitude, CellSize); //, EnvironmentalDataSource.ANY));
            if (EnviroStack.count("OceanTemp")     == 0) EnviroStack["OceanTemp"]     = EnviroData("temperature_ocean", "month", BottomLatitude, LeftmostLongitude, TopLatitude, RightmostLongitude, CellSize); //, EnvironmentalDataSource.ANY));
            if (EnviroStack.count("Precipitation") == 0) EnviroStack["Precipitation"] = EnviroData("precipitation", "month", BottomLatitude, LeftmostLongitude, TopLatitude, RightmostLongitude, CellSize); //, EnvironmentalDataSource.ANY));
            if (EnviroStack.count("FrostDays")     == 0) EnviroStack["FrostDays"]     = EnviroData("frost", "month", BottomLatitude, LeftmostLongitude, TopLatitude, RightmostLongitude, CellSize); //, EnvironmentalDataSource.ANY));
        // Initialise the cohort ID to zero
        NextCohortID = 0;
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
        //            // Create new outputs class instances (if the model is run for the whold model grid then select the grid view for the live output,
        //            // if the model is run for specific locations then use the graph view)
        //            if (SpecificLocations)
        //            {
        //
        //                // Initialise the vector of outputs instances
        //                CellOutputs = new OutputCell[_CellList.Count];
        //
        //                for (int i = 0; i < _CellList.Count; i++)
        //                {
        //                    CellOutputs[i] = new OutputCell(InitialisationFileStrings["OutputDetail"], initialisation);
        //                }
        //
        //                // Spawn a dataset viewer instance for each cell to display live model results
        //                if (initialisation.LiveOutputs)
        //                {
        //                    for (int i = 0; i < _CellList.Count; i++)
        //                    {
        //                        CellOutputs[i].SpawnDatasetViewer(NumTimeSteps);
        //                    }
        //                }
        //            
        //            }
        //            else
        //            {
        //                GridOutputs = new OutputGrid(InitialisationFileStrings["OutputDetail"], initialisation);
        //            }
        //
        //            
    }
    //----------------------------------------------------------------------------------------------
    /** \brief    Sets up the model grid within a Madingley model run
    @param initialisation An instance of the model initialisation class 
    @param scenarioParameters The parameters for the scenarios to run
    @param scenarioIndex The index of the scenario that this model is to run
     */
    void SetUpModelGrid(MadingleyModelInitialisation& initialisation) {

        //Switched order so we create cell list first then initialise cells using list rather than grid.

        unsigned NumLatCells = (unsigned) ((TopLatitude - BottomLatitude) / CellSize);
        unsigned NumLonCells = (unsigned) ((RightmostLongitude - LeftmostLongitude) / CellSize);
        //
        //                // Loop over all cells in the model
        for (unsigned ii = 0; ii < NumLatCells; ii += (unsigned) CellRarefaction) {
            for (unsigned jj = 0; jj < NumLonCells; jj += (unsigned) CellRarefaction) {
                //                        // Define a vector to hold the pair of latitude and longitude indices for this grid cell
                vector<unsigned> cellIndices(2);
                // Add the latitude and longitude indices to this vector
                cellIndices[0] = ii;
                cellIndices[1] = jj;
                // Add the vector to the list of all active grid cells
                CellList.push_back(cellIndices);
            }
        }

        EcologyTimer.Start();

        // Set up the model grid using these locations

        EcosystemModelGrid = ModelGrid(BottomLatitude, LeftmostLongitude, TopLatitude, RightmostLongitude,
                CellSize, CellSize, CellList, EnviroStack, CohortFunctionalGroupDefinitions, StockFunctionalGroupDefinitions,
                GlobalDiagnosticVariables, initialisation.TrackProcesses,  RunGridCellsInParallel);

        EcologyTimer.Stop();
        cout << "Time to initialise cells: " << EcologyTimer.GetElapsedTimeSecs() << endl;

        // Seed stocks and cohorts in the grid cells
        EcosystemModelGrid.SeedGridCellStocksAndCohorts(CellList, CohortFunctionalGroupDefinitions, StockFunctionalGroupDefinitions,
                GlobalDiagnosticVariables, NextCohortID, InitialisationFileStrings["OutputDetail"] == "high", DrawRandomly);
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
        //            if (SpecificLocations)
        //            {
        //                for (int i = 0; i < _CellList.Count; i++)
        //                {
        //                    // Set up grid cell outputs
        //                    CellOutputs[i].SetUpOutputs(EcosystemModelGrid, CohortFunctionalGroupDefinitions, StockFunctionalGroupDefinitions,
        //                        NumTimeSteps, OutputFilesSuffix, _CellList, i, EcosystemModelGrid.GetEnviroLayer("Realm", 0, _CellList[i][0], _CellList[i][1], out varExists) == 2.0);
        //
        //                    // Create initial grid cell outputs
        //                    CellOutputs[i].InitialOutputs(EcosystemModelGrid, CohortFunctionalGroupDefinitions, StockFunctionalGroupDefinitions,
        //                        _CellList, i, GlobalDiagnosticVariables, NumTimeSteps, initialisation, month, EcosystemModelGrid.GetEnviroLayer("Realm", 0, _CellList[i][0], _CellList[i][1], out varExists) == 2.0);
        //                }
        //            }
        //            else
        //            {
        //                // Set up grid outputs
        //                GridOutputs.SetupOutputs(EcosystemModelGrid, OutputFilesSuffix, NumTimeSteps, 
        //                    CohortFunctionalGroupDefinitions, StockFunctionalGroupDefinitions);
        //
        //                // Create initial grid outputs
        //                GridOutputs.InitialOutputs(EcosystemModelGrid, CohortFunctionalGroupDefinitions, StockFunctionalGroupDefinitions, _CellList, initialisation);
        //            }
        //
    }
    //----------------------------------------------------------------------------------------------
    /** \brief  Run processes for cells sequentially
     */
    void RunCellsSequentially() {
        //            // Instantiate a class to hold thread locked global diagnostic variables
        ThreadLockedParallelVariables singleThreadDiagnostics(0, 0, 0, NextCohortID);

            for (int ii = 0; ii < CellList.size(); ii++) {
                RunCell(ii, singleThreadDiagnostics,  initialisation);
            }

        // Update the variable tracking cohort unique IDs
        NextCohortID = singleThreadDiagnostics.NextCohortIDThreadLocked;
        // Take the results from the thread local variables and apply to the global diagnostic variables
        GlobalDiagnosticVariables["NumberOfCohortsExtinct"]  = singleThreadDiagnostics.Extinctions - singleThreadDiagnostics.Combinations;
        GlobalDiagnosticVariables["NumberOfCohortsProduced"] = singleThreadDiagnostics.Productions;
        GlobalDiagnosticVariables["NumberOfCohortsInModel"]  = GlobalDiagnosticVariables["NumberOfCohortsInModel"] + singleThreadDiagnostics.Productions - singleThreadDiagnostics.Extinctions;
        GlobalDiagnosticVariables["NumberOfCohortsCombined"] = singleThreadDiagnostics.Combinations;
    }
    //----------------------------------------------------------------------------------------------
    /** \brief   Run ecological processes for stocks in a specified grid cell
    @param latCellIndex The latitudinal index of the cell to run stock ecology for
    @param lonCellIndex The longitudinal index of the cell to run stock ecology for
    @param workingGridCellStocks A copy of the cohorts in the current grid cell
     * NB - acting on a copy here? needed of parallelity?
     */
    void RunWithinCellStockEcology(unsigned latCellIndex, unsigned lonCellIndex, GridCellStockHandler& workingGridCellStocks, int cellIndex) {
        // Create a local instance of the stock ecology class
        EcologyStock MadingleyEcologyStock;
        //The location of the acting stock
        vector<int> ActingStock(2);
        // Get the list of functional group indices for autotroph stocks
        vector<int> AutotrophStockFunctionalGroups = StockFunctionalGroupDefinitions.GetFunctionalGroupIndex("Heterotroph/Autotroph", "Autotroph", false);
        // Loop over autotroph functional groups
        for (int FunctionalGroup : AutotrophStockFunctionalGroups) {
            for (int ll = 0; ll < workingGridCellStocks[FunctionalGroup].size(); ll++) {
                //                    // Get the position of the acting stock
                ActingStock[0] = FunctionalGroup;
                ActingStock[1] = ll;
                // Run stock ecology
                MadingleyEcologyStock.RunWithinCellEcology(workingGridCellStocks, ActingStock, EcosystemModelGrid.GetCellEnvironment(
                        latCellIndex, lonCellIndex), EnvironmentalDataUnits, HumanNPPExtraction, StockFunctionalGroupDefinitions,
                        CurrentTimeStep, GlobalModelTimeStepUnit, ProcessTrackers[cellIndex].TrackProcesses, ProcessTrackers[cellIndex], TrackGlobalProcesses, CurrentMonth,
                        InitialisationFileStrings["OutputDetail"]);
                //workingGridCellStocks[ActingStock].TotalBiomass *= 0.75;//MB strange line - commented out in original?
            }
        }

    }
    //----------------------------------------------------------------------------------------------
    /** \brief   Run ecological processes for cohorts in a specified grid cell
    @param latCellIndex The latitudinal index of the cell to run cohort ecology for
    @param lonCellIndex The longitudinal index of the cell to run cohort ecology for
    @param workingGridCellCohorts A copy of the cohorts in the current grid cell
    @param workingGridCellStocks A copy of the stocks in the current grid cell

     * NB - acting on a copy here? needed of parallelity?
     */
    //working with *copies* of the cohorts/stocks -why? for parallellity?

    void RunWithinCellCohortEcology(unsigned latCellIndex, unsigned lonCellIndex, ThreadLockedParallelVariables& partial,
            GridCellCohortHandler& workingGridCellCohorts, GridCellStockHandler& workingGridCellStocks, string outputDetail, int cellIndex, MadingleyModelInitialisation& initialisation) {
        // Local instances of classes
        EcologyCohort MadingleyEcologyCohort;
        Activity CohortActivity;

        // A list of the original cohorts inside a particular grid cell
        vector<int> OriginalGridCellCohortsNumbers(workingGridCellCohorts.size());
        // A vector to hold the order in which cohorts will act
        vector<unsigned> RandomCohortOrder;
        // A jagged array to keep track of cohorts that are being worked on
        vector<vector<unsigned>>CohortIndices;
        // The location of the acting cohort
        vector<int> ActingCohort(2);
        // Temporary local variables
        int EcosystemModelParallelTempval1;
        int EcosystemModelParallelTempval2;
        // Boolean to pass into function to get cell environmental data to check if the specified variable exists
        bool VarExists;
        // variable to track cohort number
        unsigned TotalCohortNumber = 0;

        // Fill in the array with the number of cohorts per functional group before ecological processes are run

        for (int i = 0; i < workingGridCellCohorts.size(); i++) {
            OriginalGridCellCohortsNumbers[i] = workingGridCellCohorts[i].size();
        }

        // Initialize ecology for stocks and cohorts
        MadingleyEcologyCohort.InitializeEcology(EcosystemModelGrid.GetCellEnvironment(latCellIndex, lonCellIndex)["Cell Area"][0],
                GlobalModelTimeStepUnit, DrawRandomly);

        // Create a jagged array indexed by functional groups to hold cohort indices
        CohortIndices.resize(CohortFunctionalGroupDefinitions.GetNumberOfFunctionalGroups());

        // Loop over functional groups
        for (int ll = 0; ll < CohortFunctionalGroupDefinitions.GetNumberOfFunctionalGroups(); ll++) {
            // Dimension the number of columns in each row of the jagged array to equal number of gridCellCohorts in each functional group
            if (workingGridCellCohorts[ll].size() == 0) {
                CohortIndices[ll].push_back(0); //correct??
            } else {
                CohortIndices[ll].resize(workingGridCellCohorts[ll].size());
            }
            // Loop over gridCellCohorts in the functional group
            for (int kk = 0; kk < CohortIndices[ll].size(); kk++) {
                // Fill jagged array with indices for each cohort
                CohortIndices[ll][kk] = TotalCohortNumber;
                TotalCohortNumber += 1;
            }

        }

        if (DrawRandomly) {
            // Randomly order the cohort indices
            RandomCohortOrder = Utilities.RandomlyOrderedCohorts(TotalCohortNumber);
        } else {
            RandomCohortOrder = Utilities.NonRandomlyOrderedCohorts(TotalCohortNumber, CurrentTimeStep);
        }

        // Diagnostic biological variables don't need to be reset every cohort, but rather every grid cell
        EcosystemModelParallelTempval2 = 0;

        // Initialise eating formulations
        MadingleyEcologyCohort.EatingFormulations["Basic eating"]->InitializeEcologicalProcess(workingGridCellCohorts, workingGridCellStocks,
                CohortFunctionalGroupDefinitions, StockFunctionalGroupDefinitions, "revised predation");
        MadingleyEcologyCohort.EatingFormulations["Basic eating"]->InitializeEcologicalProcess(workingGridCellCohorts, workingGridCellStocks
                , CohortFunctionalGroupDefinitions, StockFunctionalGroupDefinitions, "revised herbivory");

        // Loop over randomly ordered gridCellCohorts to implement biological functions
        for (int ll = 0; ll < RandomCohortOrder.size(); ll++) {

            // Locate the randomly chosen cohort within the array of lists of gridCellCohorts in the grid cell
            ActingCohort = Utilities.FindJaggedArrayIndex(RandomCohortOrder[ll], CohortIndices, TotalCohortNumber);

            // Perform all biological functions except dispersal (which is cross grid cell)
            if (workingGridCellCohorts[ActingCohort[0]].size() != 0 && workingGridCellCohorts[ActingCohort].CohortAbundance > ExtinctionThreshold) {
                // Calculate number of cohorts in this functional group in this grid cell before running ecology
                EcosystemModelParallelTempval1 = workingGridCellCohorts[ActingCohort[0]].size();

                CohortActivity.AssignProportionTimeActive(workingGridCellCohorts[ActingCohort], EcosystemModelGrid.GetCellEnvironment(latCellIndex, lonCellIndex), CohortFunctionalGroupDefinitions, CurrentTimeStep, CurrentMonth);

                // Run ecology
                MadingleyEcologyCohort.RunWithinCellEcology(workingGridCellCohorts, workingGridCellStocks,
                        ActingCohort, EcosystemModelGrid.GetCellEnvironment(latCellIndex, lonCellIndex),
                        EcosystemModelGrid.GetCellDeltas(latCellIndex, lonCellIndex),
                        CohortFunctionalGroupDefinitions, StockFunctionalGroupDefinitions, CurrentTimeStep,
                        ProcessTrackers[cellIndex], partial,  outputDetail, CurrentMonth, initialisation);

                // Update the properties of the acting cohort
                MadingleyEcologyCohort.UpdateEcology(workingGridCellCohorts, workingGridCellStocks, ActingCohort,
                        EcosystemModelGrid.GetCellEnvironment(latCellIndex, lonCellIndex), EcosystemModelGrid.GetCellDeltas(
                        latCellIndex, lonCellIndex), CohortFunctionalGroupDefinitions, StockFunctionalGroupDefinitions, CurrentTimeStep,
                        ProcessTrackers[cellIndex]);

                // Add newly produced cohorts to the tracking variable
                EcosystemModelParallelTempval2 += workingGridCellCohorts[ActingCohort[0]].size() - EcosystemModelParallelTempval1;


                // Check that the mass of individuals in this cohort is still >= 0 after running ecology
                assert(workingGridCellCohorts[ActingCohort].IndividualBodyMass >= 0.0 && "Biomass < 0 for this cohort");
            }

            // Check that the mass of individuals in this cohort is still >= 0 after running ecology
            if (workingGridCellCohorts[ActingCohort[0]].size() > 0)assert(workingGridCellCohorts[ActingCohort].IndividualBodyMass >= 0.0 && "Biomass < 0 for this cohort");
        }


        // Update diagnostics of productions
        partial.Productions += EcosystemModelParallelTempval2;

        RunExtinction(latCellIndex, lonCellIndex, partial, workingGridCellCohorts, cellIndex);



        // Merge cohorts, if necessary
        if (workingGridCellCohorts.GetNumberOfCohorts() > initialisation.MaxNumberOfCohorts) {
            partial.Combinations = CohortMerger.MergeToReachThresholdFast(workingGridCellCohorts, workingGridCellCohorts.GetNumberOfCohorts(), initialisation.MaxNumberOfCohorts);

            //Run extinction a second time to remove those cohorts that have been set to zero abundance when merging
            RunExtinction(latCellIndex, lonCellIndex, partial, workingGridCellCohorts, cellIndex);
        } else
            partial.Combinations = 0;

        // Write out the updated cohort numbers after all ecological processes have occured
        EcosystemModelGrid.SetGridCellCohorts(workingGridCellCohorts, latCellIndex, lonCellIndex);
    }

    //----------------------------------------------------------------------------------------------
    /** \brief Carries out extinction on cohorts that have an abundance below a defined extinction threshold */
    void RunExtinction(unsigned latCellIndex, unsigned lonCellIndex, ThreadLockedParallelVariables& partial,
            GridCellCohortHandler workingGridCellCohorts, int cellIndex) {
        bool VarExists;

        // Loop over cohorts and remove any whose abundance is below the extinction threshold
        for (int kk = 0; kk < CohortFunctionalGroupDefinitions.GetNumberOfFunctionalGroups(); kk++) {
            // Create a list to hold the cohorts to remove
            vector<int> CohortIndicesToRemove;

            // Loop through each cohort in the functional group
            for (int ll = 0; ll < workingGridCellCohorts[kk].size(); ll++) {
                // If cohort abundance is less than the extinction threshold then add to the list for extinction
                if (workingGridCellCohorts[kk][ll].CohortAbundance > ExtinctionThreshold || workingGridCellCohorts[kk][ll].IndividualBodyMass <= 1.e-300) {
                    CohortIndicesToRemove.push_back(ll);

                    partial.Extinctions += 1;

                    //                        // If track processes is set and output detail is set to high and the cohort being made extinct has never been merged,
                    //                        // then output its mortality profile
                    //                        if (ProcessTrackers[cellIndex].TrackProcesses && (InitialisationFileStrings["OutputDetail"] == "high") && (workingGridCellCohorts[kk][ll].CohortID.Count == 1))
                    //                        {
                    //                            ProcessTrackers[cellIndex].OutputMortalityProfile(workingGridCellCohorts[kk][ll].CohortID[0]);
                    //                        }
                }
            }

            // Code to add the biomass to the biomass pool and dispose of the cohort
            for (int ll = (CohortIndicesToRemove.size() - 1); ll >= 0; ll--) {
                //                    // Add biomass of the extinct cohort to the organic matter pool
                EcosystemModelGrid.SetEnviroLayer("Organic Pool", 0, EcosystemModelGrid.GetEnviroLayer("Organic Pool", 0, latCellIndex, lonCellIndex, VarExists) +
                        (workingGridCellCohorts[kk][CohortIndicesToRemove[ll]].IndividualBodyMass + workingGridCellCohorts[kk][CohortIndicesToRemove[ll]].IndividualReproductivePotentialMass) * workingGridCellCohorts[kk][CohortIndicesToRemove[ll]].CohortAbundance, latCellIndex, lonCellIndex);
                assert(EcosystemModelGrid.GetEnviroLayer("Organic Pool", 0, latCellIndex, lonCellIndex, VarExists) >= 0 && "Organic pool < 0");
                //
                //                    if (ProcessTrackers[cellIndex].TrackProcesses && SpecificLocations == true)
                //                        ProcessTrackers[cellIndex].RecordExtinction(latCellIndex, lonCellIndex, CurrentTimeStep, workingGridCellCohorts[kk][CohortIndicesToRemove[ll]].Merged, workingGridCellCohorts[kk][CohortIndicesToRemove[ll]].CohortID);

                // Remove the extinct cohort from the list of cohorts
                workingGridCellCohorts[kk].erase(workingGridCellCohorts[kk].begin() + CohortIndicesToRemove[ll]);


            }

        }

    }
    //----------------------------------------------------------------------------------------------
    /** \brief Run ecological processes that operate across grid cells */
    void RunCrossGridCellEcology(unsigned& dispersals) {
        // Loop through each grid cell, and run dispersal for each.
        // Note that currently dispersal is not parallelised, although it could be (though care would need to be taken to ensure that necessary variables are thread-locked
        
        for (int ii = 0; ii < CellList.size(); ii++) {

            // We have looped through individual cells and calculated ecological processes for each. Now do this for cross grid cell dispersal

            disperser.RunCrossGridCellEcologicalProcess(CellList[ii], EcosystemModelGrid,  CohortFunctionalGroupDefinitions, StockFunctionalGroupDefinitions, CurrentMonth);

        }
        // Apply the changes from dispersal
        disperser.UpdateCrossGridCellEcology(EcosystemModelGrid, dispersals);

    }
    //----------------------------------------------------------------------------------------------
    /** \brief Make a record of the properties of the intial model cohorts in the new cohorts output file */
    void RecordInitialCohorts() {
        //            int i = 0;
        //            foreach (uint[] cell in _CellList)
        //            {
        //                if (ProcessTrackers[i].TrackProcesses)
        //                {
        //
        //                    GridCellCohortHandler TempCohorts = EcosystemModelGrid.GetGridCellCohorts(cell[0], cell[1]);
        //
        //                    for (int FunctionalGroup = 0; FunctionalGroup < TempCohorts.Count; FunctionalGroup++)
        //                    {
        //                        foreach (Cohort item in TempCohorts[FunctionalGroup])
        //                        {
        //                            ProcessTrackers[i].RecordNewCohort(cell[0], cell[1], 0, item.CohortAbundance, item.AdultMass, item.FunctionalGroupIndex,
        //                                new List<uint> { uint.MaxValue }, item.CohortID[0]);
        //                        }
        //                    }
        //                }
        //                i += 1;
        //            }
        //        }
        //
    }
    //----------------------------------------------------------------------------------------------

};
#endif