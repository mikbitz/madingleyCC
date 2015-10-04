#ifndef MADINGLEYMODEL_H
#define MADINGLEYMODEL_H
#include <vector>
#include <map>
#include <ScenarioParameterInitialisation.h>
#include <MadingleyModelInitialisation.h>
#include <FunctionalGroupDefinitions.h>
#include <Stopwatch.h>
#include <CohortMerge.h>
#include <ModelGrid.h>
#include <EcologyCrossGridCell.h>
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
//Changes
//MB 7/8/2015
//change uint to unsigned
//change _CellList to vector <vector <unsigned> >
//change sortedLists to map<>
//namespace Madingley
//{   
//
/** \brief Thread-local variables for tracking extinction and production of cohorts

//    /// <todo>Needs a little tidying and checking of access levels</todo>
*/

//
/** \brief The ecosystem model */
class MadingleyModel
    {
    public:
    private:
//        
/** An instance of the cohort functional group definitions for this model */
     FunctionalGroupDefinitions CohortFunctionalGroupDefinitions;
/** An instance of the stock functional group definitions for this model */
     FunctionalGroupDefinitions StockFunctionalGroupDefinitions;
//
/** \brief
 A list of environmental data layers
*/
     map<string, EnviroData> EnviroStack;
//
/** \brief An instance of ModelGrid to hold the grid to be used in this model */
     ModelGrid EcosystemModelGrid;
/** \brief An instance of the cross grid cell ecology class */
     EcologyCrossGridCell MadingleyEcologyCrossGridCell;

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

/** \brief Whether the model will be run for specific locations, instead of for the whole model grid */
        bool SpecificLocations;        
/** \brief An instance of StopWatch to time individual time steps */
        StopWatch TimeStepTimer;
        StopWatch EcologyTimer;
        StopWatch OutputTimer;
/** \brief
 An array of instances of the output class to deal with grid cell outputs
*/
//        private OutputCell[] CellOutputs;
//
/** \brief
 An array of indices of process trackers for each grid cell
*/
vector<ProcessTracker> ProcessTrackers;
//
/** \brief
 An instance of a global process tracker to track across the model grid
*/
GlobalProcessTracker TrackGlobalProcesses;
//
/** \brief An instance of the output class to deal with global outputs */
//        private OutputGlobal GlobalOutputs;
/** \brief
 An instance of the output class to deal with gridded outputs
*/
//        private OutputGrid GridOutputs;
//
/** \brief The suffix to be applied to files output by this model instance */
        string OutputFilesSuffix;
/** \brief
 A sorted list of strings from the initialisation file
*/
          map<string, string> InitialisationFileStrings;
/** \brief A sorted list of strings for environmental data units */
          map<string, string> EnvironmentalDataUnits;
/** \brief Get the human appropriation of NPP scenario to use */
          string HumanNPPExtraction;
/** A variable to increment for the purposes of giving each cohort a unique ID */
          long long NextCohortID;
//     
/** \brief Variable to track the number of cohorts that have dispersed. Doesn't need to be thread-local because all threads have converged prior to running cross-grid-cell processes */
         unsigned Dispersals;
//
/** \brief Instance of the class to perform general functions */
          UtilityFunctions Utilities;
//
/** \brief An instance of the merging class */
          CohortMerge CohortMerger;
/** \brief
 Initializes the ecosystem model

@param initialisation An instance of the model initialisation class 
@param scenarioParameters The parameters for the scenarios to run
@param scenarioIndex The index of the scenario that this model is to run
@param outputFilesSuffix The suffix to be applied to all outputs from this model run
@param globalModelTimeStepUnit The time step unit used in the model
*/
    public:
     MadingleyModel(MadingleyModelInitialisation& initialisation, ScenarioParameterInitialisation& scenarioParameters, int scenarioIndex,
            string outputFilesSuffix, string globalModelTimeStepUnit)
        {         
//            // Assign the properties for this model run
            AssignModelRunProperties(initialisation, scenarioParameters, scenarioIndex, outputFilesSuffix);
//
//            // Set up list of global diagnostics
            SetUpGlobalDiagnosticsList();

//Class level variable - currently doesn't need initializing here
//            // Initialise the cell list
//            _CellList = new List<uint[]>();
//
//            // Set up the model grid

            SetUpModelGrid(initialisation, scenarioParameters, scenarioIndex);
//
//            // Set up model outputs
            SetUpOutputs(initialisation);
//
//            // Make the initial outputs
            InitialOutputs(outputFilesSuffix, initialisation, CurrentMonth);
//
//            // Instance the array of process trackers
//            ProcessTrackers = new ProcessTracker[_CellList.Count];
//
//            // Temporary variables
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
//            if (SpecificLocations) initialisation.RunRealm = "";
//
//            // Record the initial cohorts in the process trackers
            RecordInitialCohorts();
//Currently don't need to initalise these as they are created at class level
//            // Initialise the class for cross-grid-cell ecology 
           // MadingleyEcologyCrossGridCell = new EcologyCrossGridCell();
//
//            // Initialise the time step timer
            //TimeStepTimer = new StopWatch();
            //EcologyTimer = new StopWatch();
            //OutputTimer = new StopWatch();
//
            // Set the global model time step unit
            GlobalModelTimeStepUnit = globalModelTimeStepUnit;

//            // Initialise the utility functions
//            Utilities = new UtilityFunctions();
//end of initialisations
            
//            // Initialise the cohort merger - this is just to set where the random seed comes from
            CohortMerger.SetRandom(DrawRandomly);
//
        }
//
/** \brief
 Run the global ecosystem model

@param initialisation The initialization details for the current set of model simulations
*/
void RunMadingley(MadingleyModelInitialisation& initialisation)
        {            
            // Write out model run details to the console
            cout<<"Running model"<<endl;
            cout<<"Number of time steps is: "<<NumTimeSteps<<endl;

//            // Temporary variable
            bool varExists;
            Dispersals=0;
//            
//
//             // Run the model
             for (unsigned hh = 0; hh < NumTimeSteps; hh += 1)
             {
                 cout<<"Running time step "<<hh+1<<"..."<<endl;
//
//                 // Start the timer
                 TimeStepTimer.Start();

                 // Get current time step and month
                 CurrentTimeStep = hh;
                 CurrentMonth = Utilities.GetCurrentMonth(hh,GlobalModelTimeStepUnit);
//
//                 // Initialise cross grid cell ecology
                 MadingleyEcologyCrossGridCell.InitializeCrossGridCellEcology(GlobalModelTimeStepUnit, DrawRandomly, initialisation);
//
                 EcologyTimer.Start();
//
//                 // Loop over grid cells and run biological processes
//                 if (RunGridCellsInParallel)
//                 {
//                     // Run cells in parallel
//                     RunCellsInParallel(initialisation);
//                 }
//                 else
//                 {
//                     // Run cells in sequence
                     RunCellsSequentially(initialisation);
//                 }
//
//                 
//
                 EcologyTimer.Stop();
                 cout<<"Within grid ecology took: "<<EcologyTimer.GetElapsedTimeSecs()<<endl;
//                 // Run the garbage collector. Note that it works in the background so may take a little while
//                 // Needs to be done to ensure cohorts are deleted properly
//                 GC.Collect();
//
//                 if(TrackGlobalProcesses.TrackProcesses) TrackGlobalProcesses.StoreNPPGrid(hh);
//
                 EcologyTimer.Start();
//                 // Run cross grid cell ecology
                 RunCrossGridCellEcology(Dispersals, initialisation.DispersalOnly, initialisation);
//
                 EcologyTimer.Stop();
                 cout<<"Across grid ecology took: "<< EcologyTimer.GetElapsedTimeSecs()<<endl;
//
//                 // Run the garbage collector. Note that it works in the background so may take a little while
//                 // Needs to be done here to ensure cohorts are deleted properly
//                 GC.Collect();
//                 
//                 // Stop the timer
                 TimeStepTimer.Stop();
//
                 OutputTimer.Start();
//
//                 // Write the global outputs for this time step
//                 GlobalOutputs.TimeStepOutputs(EcosystemModelGrid, CurrentTimeStep, CurrentMonth, TimeStepTimer,CohortFunctionalGroupDefinitions,
//                     StockFunctionalGroupDefinitions,_CellList,GlobalDiagnosticVariables, initialisation);
//
                 OutputTimer.Stop();
                 cout<<"Global Outputs took: "<<OutputTimer.GetElapsedTimeSecs()<<endl;
//
//
                 OutputTimer.Start();
//
//                 if (SpecificLocations)
//                 {
//                     // Loop over grid cells and write outputs
//                     for (int i = 0; i < _CellList.Count; i++).
//                ToArray();
//                     {
//                         // Write out the grid cell outputs for this time step
//                         CellOutputs[i].TimeStepOutputs(EcosystemModelGrid, CohortFunctionalGroupDefinitions, StockFunctionalGroupDefinitions,
//                             _CellList, i, GlobalDiagnosticVariables, TimeStepTimer, NumTimeSteps, CurrentTimeStep, initialisation, CurrentMonth, EcosystemModelGrid.GetEnviroLayer("Realm", 0, _CellList[i][0], _CellList[i][1], out varExists) == 2.0);
//                         
//                         // Write out trophic flow data for this time step
//                         if(ProcessTrackers[i].TrackProcesses) ProcessTrackers[i].WriteTimeStepTrophicFlows(CurrentTimeStep, EcosystemModelGrid.NumLatCells, EcosystemModelGrid.NumLonCells, initialisation,
//                             EcosystemModelGrid.GetEnviroLayer("Realm", 0, _CellList[i][0], _CellList[i][1], out varExists) == 2.0);
//                 
//                     }
//                 }
//                 else
//                 {
//                     // Write out grid outputs for this time step
//                     GridOutputs.TimeStepOutputs(EcosystemModelGrid, CohortFunctionalGroupDefinitions, StockFunctionalGroupDefinitions, _CellList,
//                         CurrentTimeStep, initialisation);
//                 }
//
//
                 OutputTimer.Stop();
                 cout<<"Cell/Grid Outputs took: "<<OutputTimer.GetElapsedTimeSecs()<<endl;

                 // Write the results of dispersal to the console
                 cout<<"Number of cohorts that dispersed this time step: "<<Dispersals<<endl;

                 //Dispersals = 0;
//   
//                 
//
//
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
//
//            if (SpecificLocations)
//            {
//                // Loop over grid cells and write the final grid cell outputs
//                for (int i = 0; i < _CellList.Count; i++)
//                {
//                    CellOutputs[i].FinalOutputs(EcosystemModelGrid, CohortFunctionalGroupDefinitions, StockFunctionalGroupDefinitions,
//                        _CellList, i, GlobalDiagnosticVariables, initialisation, CurrentMonth, EcosystemModelGrid.GetEnviroLayer("Realm", 0, _CellList[i][0], _CellList[i][1], out varExists) == 2.0);
//                }
//            }
//            else
//            {
//                // Write the final grid outputs
//                GridOutputs.FinalOutputs();
//            }
//
//            
        }
//
/** \brief
 A method to run the main ecosystem model loop in parallel (latitudinal strips)

@param latCellIndex The latitude index of the cell to run
@param lonCellIndex The longitude index of the cell to run
@param partial A threadlockedparallelvariable that is used to pass global diagnostic information back with locking or race conditions
@param dispersalOnly Whether to run dispersal only (i.e. to turn all other ecological processes off
 <remarks>Note that variables and instances of classes that are written to within this method MUST be local within this method to prevent 
 race issues and multiple threads attempting to write to the same variable when running the program in parallel</remarks>
*/
void RunCell(int cellIndex, ThreadLockedParallelVariables& partial, bool dispersalOnly, MadingleyModelInitialisation& initialisation)
        {
//
//            // Create a temporary internal copy of the grid cell cohorts
            GridCellCohortHandler WorkingGridCellCohorts = EcosystemModelGrid.GetGridCellCohorts(CellList[cellIndex][0], CellList[cellIndex][1]);
//
//            // Create a temporary internal copy of the grid cell stocks
            GridCellStockHandler WorkingGridCellStocks = EcosystemModelGrid.GetGridCellStocks(CellList[cellIndex][0], CellList[cellIndex][1]);
//
//            // Run stock ecology
            RunWithinCellStockEcology(CellList[cellIndex][0], CellList[cellIndex][1], WorkingGridCellStocks, cellIndex);
//
//            // Run within cell ecology if we are not doing dispersal only
//            if (dispersalOnly)
//            {
//                // Run cohort ecology
//                RunWithinCellDispersalOnly(_CellList[cellIndex][0], _CellList[cellIndex][1], partial, WorkingGridCellCohorts, WorkingGridCellStocks);
//            }
//            else
//            {
//                // Run cohort ecology
                RunWithinCellCohortEcology(CellList[cellIndex][0],CellList[cellIndex][1], partial, WorkingGridCellCohorts, WorkingGridCellStocks, InitialisationFileStrings["OutputDetail"], cellIndex, initialisation);
//
//            }
//
//            // For runs with specific locations and where track processes has been specified, write out mass flows data and reset the mass flow tracker 
//            // for the next time step
//            if (SpecificLocations && ProcessTrackers[cellIndex].TrackProcesses)
//            {
//                ProcessTrackers[cellIndex].EndTimeStepPredationTracking(CurrentTimeStep);
//                ProcessTrackers[cellIndex].EndTimeStepHerbvioryTracking(CurrentTimeStep);
//            }
//
        }
//
/** \brief
 Read in the specified locations in which to run the model

@param specificLocationsFile The name of the file with specific locations information
@param outputPath The path to the output folder in which to copy the specific locations file
*/
void ReadSpecificLocations(string specificLocationsFile, string outputPath)
        {
//
//            List<double> LatitudeList = new List<double>();
//            List<double> LongitudeList = new List<double>();
//
//            Console.WriteLine("Reading in specific location data");
//            Console.WriteLine("");
//
//            // construct file name
//            string FileString = "msds:csv?file=input/Model setup/" + specificLocationsFile + "&openMode=readOnly";
//
//            // Read in the data
//            DataSet InternalData = DataSet.Open(FileString);
//
//            foreach (Variable v in InternalData.Variables)
//            {
//                //Get the name of the variable currently referenced in the dataset
//                string HeaderName = v.Name;
//                //Copy the values for this variable into an array
//                var TempValues = v.GetData();
//
//                switch (HeaderName.ToLower())
//                {
//                    // Add the latitude and longitude values to the appropriate list
//                    case "latitude":
//                        for (int ii = 0; ii < TempValues.Length; ii++) LatitudeList.Add(Convert.ToDouble(TempValues.GetValue(ii).ToString()));
//                        break;
//                    case "longitude":
//                        for (int ii = 0; ii < TempValues.Length; ii++) LongitudeList.Add(Convert.ToDouble(TempValues.GetValue(ii).ToString()));
//                        break;
//                    default:
//                        Console.WriteLine("Variable defined in the specific location file but not processed: ", HeaderName);
//                        break;
//                }
//            }
//
//            // Loop over cells defined in the specific locations file
//            for (int ii = 0; ii < LatitudeList.Count; ii++)
//            {
//                // Define a vector to hold the longitude and latitude index for this cell
//                uint[] cellIndices = new uint[2];
//
//                // Get the longitude and latitude indices for the current grid cell
//                cellIndices[0] = (uint)Math.Floor((LatitudeList.ElementAt(ii) - BottomLatitude) / CellSize);
//                cellIndices[1] = (uint)Math.Floor((LongitudeList.ElementAt(ii) - LeftmostLongitude) / CellSize);
//
//                // Add these indices to the list of active cells
//                _CellList.Add(cellIndices);
//            }
//            
//
        }

/** \brief
 Assigns the properties of the current model run

@param initialisation An instance of the model initialisation class 
@param scenarioParameters The parameters for the scenarios to run
@param scenarioIndex The index of the scenario that this model is to run
@param outputFilesSuffix The suffix to be applied to all outputs from this model run
*/
void AssignModelRunProperties(MadingleyModelInitialisation& initialisation, ScenarioParameterInitialisation& scenarioParameters, int scenarioIndex,
            string outputFilesSuffix)
        {
//            // Assign the properties of this model run from the same properties in the specified model initialisation
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
             HumanNPPExtraction = scenarioParameters.scenarioParameters["Human NPP Extraction"][scenarioIndex];
             OutputFilesSuffix = outputFilesSuffix;
             EnvironmentalDataUnits = initialisation.Units;
//
//            // Set whether the model run is for specific locations
//            if (InitialisationFileStrings.ContainsKey("Locations")) SpecificLocations = true;
//
//            // If the model run is for a whole grid, then create envirodata classes with the environmental data
            if (!SpecificLocations){
                // TEMP COMMENT OUT OF TEMPERATURE
                //MB originally these were done with "FetchClimate" but this needs to be replaced
                if (EnviroStack.count("Temperature")  ==0) EnviroStack["Temperature"]  =EnviroData("temperature", "month", BottomLatitude, LeftmostLongitude, TopLatitude, RightmostLongitude, CellSize);//, EnvironmentalDataSource.ANY));
                if (EnviroStack.count("LandDTR")      ==0) EnviroStack["LandDTR"]      =EnviroData("land_dtr", "month", BottomLatitude, LeftmostLongitude, TopLatitude, RightmostLongitude, CellSize);//, EnvironmentalDataSource.ANY));
                if (EnviroStack.count("OceanTemp")    ==0) EnviroStack["OceanTemp"]    =EnviroData("temperature_ocean", "month", BottomLatitude, LeftmostLongitude, TopLatitude, RightmostLongitude, CellSize);//, EnvironmentalDataSource.ANY));
                if (EnviroStack.count("Precipitation")==0) EnviroStack["Precipitation"]=EnviroData("precipitation", "month", BottomLatitude, LeftmostLongitude, TopLatitude, RightmostLongitude, CellSize);//, EnvironmentalDataSource.ANY));
                if (EnviroStack.count("FrostDays")    ==0) EnviroStack["FrostDays"]    =EnviroData("frost", "month", BottomLatitude, LeftmostLongitude, TopLatitude, RightmostLongitude, CellSize);//, EnvironmentalDataSource.ANY));
            }
//
//            // Initialise the cohort ID to zero
            NextCohortID = 0;
//
        }

/** \brief
 Sets up the list of global diagnostic variables
*/
void SetUpGlobalDiagnosticsList()
        {
//            //Instantiate the global diagnostic variables
//            GlobalDiagnosticVariables = new SortedList<string, double>();
            // Add global diagnostic variables
            GlobalDiagnosticVariables["NumberOfCohortsExtinct"]= 0.0;
            GlobalDiagnosticVariables["NumberOfCohortsProduced"]= 0.0;
            GlobalDiagnosticVariables["NumberOfCohortsCombined"]= 0.0;
            GlobalDiagnosticVariables["NumberOfCohortsInModel"]= 0.0;
            GlobalDiagnosticVariables["NumberOfStocksInModel"]= 0.0;

        }

/** \brief
 Sets up the model outputs

@param initialisation An instance of the model initialisation class
*/
void SetUpOutputs(MadingleyModelInitialisation initialisation)
        {
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
//
//                // Spawn dataset viewer to display live grid results
//                if (initialisation.LiveOutputs)
//                {
//                    GridOutputs.SpawnDatasetViewer();
//                }
//            }
//
//            
        }

/** \brief
 Sets up the model grid within a Madingley model run

@param initialisation An instance of the model initialisation class 
@param scenarioParameters The parameters for the scenarios to run
@param scenarioIndex The index of the scenario that this model is to run
*/
void SetUpModelGrid(MadingleyModelInitialisation& initialisation, ScenarioParameterInitialisation& scenarioParameters, int scenarioIndex)
        {
//            // If the intialisation file contains a column pointing to another file of specific locations, and if this column is not blank then read the 
//            // file indicated
//            if (SpecificLocations)
//            {
//                // Read the file containing list of specific locations to run the model for (this method also assigns the latitude and longitude indices)
//                if (!InitialisationFileStrings["Locations"].Equals("")) this.ReadSpecificLocations(InitialisationFileStrings["Locations"], initialisation.OutputPath);
//                // Set up the model grid using these locations
//                EcosystemModelGrid = new ModelGrid(BottomLatitude, LeftmostLongitude, TopLatitude, RightmostLongitude,
//                    CellSize, CellSize, _CellList, EnviroStack, CohortFunctionalGroupDefinitions, StockFunctionalGroupDefinitions,
//                    GlobalDiagnosticVariables, initialisation.TrackProcesses, SpecificLocations,RunGridCellsInParallel);
//
//            }
//            else
//            {
//
//                //Switched order so we create cell list first then initialise cells using list rather than grid.
//
                unsigned NumLatCells = (unsigned)((TopLatitude - BottomLatitude) / CellSize);
                unsigned NumLonCells = (unsigned)((RightmostLongitude - LeftmostLongitude) / CellSize);
//
//                // Loop over all cells in the model
                for (unsigned ii = 0; ii < NumLatCells; ii += (unsigned)CellRarefaction)
                {
                    for (unsigned jj = 0; jj < NumLonCells; jj += (unsigned)CellRarefaction)
                    {
//                        // Define a vector to hold the pair of latitude and longitude indices for this grid cell
                        vector<unsigned> cellIndices(2);
//
//                        // Add the latitude and longitude indices to this vector
                        cellIndices[0] = ii;
                        cellIndices[1] = jj;
//
//                        // Add the vector to the list of all active grid cells
                        CellList.push_back(cellIndices);
//
                    }
                }
//
//                EcologyTimer = new StopWatch();
                EcologyTimer.Start();

//
//                // Set up a full model grid (i.e. not for specific locations)
//                // COMMENTING OUT OLD CODE FOR INITIALISING A MODEL GRID
//                //EcosystemModelGrid = new ModelGrid(BottomLatitude, LeftmostLongitude, TopLatitude, RightmostLongitude,
//                //    CellSize, CellSize, CellRarefaction, EnviroStack, CohortFunctionalGroupDefinitions, StockFunctionalGroupDefinitions, 
//                //    GlobalDiagnosticVariables, ref NextCohortID, initialisation.TrackProcesses,DrawRandomly, SpecificLocations);
//
//                // Set up the model grid using these locations
                
                EcosystemModelGrid = ModelGrid(BottomLatitude, LeftmostLongitude, TopLatitude, RightmostLongitude,
                    CellSize, CellSize, CellList, EnviroStack, CohortFunctionalGroupDefinitions, StockFunctionalGroupDefinitions,
                    GlobalDiagnosticVariables, initialisation.TrackProcesses, SpecificLocations, RunGridCellsInParallel);

                EcologyTimer.Stop();
              cout<<"Time to initialise cells: "<< EcologyTimer.GetElapsedTimeSecs()<<endl;

//
//                Console.ForegroundColor = ConsoleColor.Red;
//                Console.WriteLine("Madingley Model memory usage post grid cell seed: {0}", GC.GetTotalMemory(true) / 1E9, " (G Bytes)\n");
//                Console.ForegroundColor = ConsoleColor.White;
//
//            }
//
//            // When the last simulation for the current scenario
//            if ((scenarioParameters.scenarioSimulationsNumber.Count == 1) && (scenarioIndex == scenarioParameters.scenarioSimulationsNumber[scenarioIndex] - 1)) EnviroStack.Clear();
//
//            // Seed stocks and cohorts in the grid cells
            EcosystemModelGrid.SeedGridCellStocksAndCohorts(CellList, CohortFunctionalGroupDefinitions, StockFunctionalGroupDefinitions, 
                GlobalDiagnosticVariables, NextCohortID, InitialisationFileStrings["OutputDetail"]=="high",DrawRandomly, initialisation.DispersalOnly, InitialisationFileStrings["DispersalOnlyType"]); 
//
//            Console.ForegroundColor = ConsoleColor.Red;
//            Console.WriteLine("Madingley Model memory usage pre Collect: {0}", Math.Round(GC.GetTotalMemory(true) / 1E9, 2), " (GBytes)");
//            Console.ForegroundColor = ConsoleColor.White;
//            GC.Collect();
//
//            Console.ForegroundColor = ConsoleColor.Red;
//            Console.WriteLine("Madingley Model memory usage post Collect: {0}", Math.Round(GC.GetTotalMemory(true) / 1E9, 5), " (GBytes)\n");
//            Console.ForegroundColor = ConsoleColor.White;
//
        }

/** \brief
 Generates the initial outputs for this model run

@param outputFilesSuffix The suffix to be applied to all outputs from this model run
*/
void InitialOutputs(string outputFilesSuffix, MadingleyModelInitialisation initialisation, unsigned month)
        {
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

/** \brief Run processes for cells in parallel */
//        public void RunCellsInParallel(MadingleyModelInitialisation initialisation)
//        {
//
//            // Create temporary variables to hold extinctions and productions in the parallel loop;
//            int extinctions = 0, productions = 0, combinations = 0;
//
//            if (initialisation.RunRealm == "")
//            {
//                // Run a parallel loop over rows
//                Parallel.For(0, _CellList.Count, () => new ThreadLockedParallelVariables { Extinctions = 0, Productions = 0, Combinations = 0, NextCohortIDThreadLocked = NextCohortID }, (ii, loop, threadTrackedDiagnostics) =>
//                {
//                    RunCell(ii, threadTrackedDiagnostics, initialisation.DispersalOnly, initialisation);
//                    return threadTrackedDiagnostics;
//                },
//                 (threadTrackedDiagnostics) =>
//                 {
//                     Interlocked.Add(ref extinctions, threadTrackedDiagnostics.Extinctions);
//                     Interlocked.Add(ref productions, threadTrackedDiagnostics.Productions);
//                     Interlocked.Add(ref combinations, threadTrackedDiagnostics.Combinations);
//                     Interlocked.Exchange(ref NextCohortID, threadTrackedDiagnostics.NextCohortIDThreadLocked);
//                 }
//                 );
//
//
//            }
//            else
//            {
//
//                if (initialisation.RunRealm == "marine")
//                {
//
//                    // Run a parallel loop over rows
//                    Parallel.For(0, _CellList.Count, () => new ThreadLockedParallelVariables { Extinctions = 0, Productions = 0, Combinations = 0, NextCohortIDThreadLocked = NextCohortID }, (ii, loop, threadTrackedDiagnostics) =>
//                    {
//                        if (EcosystemModelGrid.GetCellEnvironment(_CellList[ii][0], _CellList[ii][1])["Realm"][0] == 2.0) RunCell(ii, threadTrackedDiagnostics, initialisation.DispersalOnly, initialisation);
//                        return threadTrackedDiagnostics;
//                    },
//                     (threadTrackedDiagnostics) =>
//                     {
//                         Interlocked.Add(ref extinctions, threadTrackedDiagnostics.Extinctions);
//                         Interlocked.Add(ref productions, threadTrackedDiagnostics.Productions);
//                         Interlocked.Add(ref combinations, threadTrackedDiagnostics.Combinations);
//                         Interlocked.Exchange(ref NextCohortID, threadTrackedDiagnostics.NextCohortIDThreadLocked);
//                     }
//                     );
//                }
//                else
//                {
//
//                    // Run a parallel loop over rows
//                    Parallel.For(0, _CellList.Count, () => new ThreadLockedParallelVariables { Extinctions = 0, Productions = 0, Combinations = 0, NextCohortIDThreadLocked = NextCohortID }, (ii, loop, threadTrackedDiagnostics) =>
//                    {
//                        if (EcosystemModelGrid.GetCellEnvironment(_CellList[ii][0], _CellList[ii][1])["Realm"][0] == 1.0) RunCell(ii, threadTrackedDiagnostics, initialisation.DispersalOnly, initialisation);
//                        return threadTrackedDiagnostics;
//                    },
//                     (threadTrackedDiagnostics) =>
//                     {
//                         Interlocked.Add(ref extinctions, threadTrackedDiagnostics.Extinctions);
//                         Interlocked.Add(ref productions, threadTrackedDiagnostics.Productions);
//                         Interlocked.Add(ref combinations, threadTrackedDiagnostics.Combinations);
//                         Interlocked.Exchange(ref NextCohortID, threadTrackedDiagnostics.NextCohortIDThreadLocked);
//                     }
//                     );
//                }
//            }
//
//            Console.WriteLine("\n");
//
//            // Take the results from the thread local variables and apply to the global diagnostic variables
//            GlobalDiagnosticVariables["NumberOfCohortsExtinct"] = extinctions - combinations;
//            GlobalDiagnosticVariables["NumberOfCohortsProduced"] = productions;
//            GlobalDiagnosticVariables["NumberOfCohortsInModel"] = GlobalDiagnosticVariables["NumberOfCohortsInModel"] + productions - extinctions;
//            GlobalDiagnosticVariables["NumberOfCohortsCombined"] = combinations;
//        }
//
/** \brief
 Run processes for cells sequentially
*/
void RunCellsSequentially(MadingleyModelInitialisation& initialisation)
        {
//            // Instantiate a class to hold thread locked global diagnostic variables
            ThreadLockedParallelVariables singleThreadDiagnostics(0,0,0,NextCohortID);

            if (initialisation.RunRealm == "")
            {

                for (int ii = 0; ii < CellList.size(); ii++)
                {
                    RunCell(ii, singleThreadDiagnostics, initialisation.DispersalOnly, initialisation);
                }
            }

            else
            {

                if (initialisation.RunRealm == "marine")
                {
                    for (int ii = 0; ii < CellList.size(); ii++)
                    {
                        if (EcosystemModelGrid.GetCellEnvironment(CellList[ii][0], CellList[ii][1])["Realm"][0] == 2.0) RunCell(ii, singleThreadDiagnostics, initialisation.DispersalOnly, initialisation);
                    }
                }
                else
                {
                    for (int ii = 0; ii < CellList.size(); ii++)
                    {
                        if (EcosystemModelGrid.GetCellEnvironment(CellList[ii][0], CellList[ii][1])["Realm"][0] == 1.0) RunCell(ii, singleThreadDiagnostics, initialisation.DispersalOnly, initialisation);
                    }
                }
            }

//            // Update the variable tracking cohort unique IDs
            NextCohortID = singleThreadDiagnostics.NextCohortIDThreadLocked;
//
//            // Take the results from the thread local variables and apply to the global diagnostic variables
            GlobalDiagnosticVariables["NumberOfCohortsExtinct"] = singleThreadDiagnostics.Extinctions - singleThreadDiagnostics.Combinations;
            GlobalDiagnosticVariables["NumberOfCohortsProduced"] = singleThreadDiagnostics.Productions;
            GlobalDiagnosticVariables["NumberOfCohortsInModel"] = GlobalDiagnosticVariables["NumberOfCohortsInModel"] + singleThreadDiagnostics.Productions - singleThreadDiagnostics.Extinctions;
            GlobalDiagnosticVariables["NumberOfCohortsCombined"] = singleThreadDiagnostics.Combinations;
        }

/** \brief
 Run ecological processes for stocks in a specified grid cell

@param latCellIndex The latitudinal index of the cell to run stock ecology for
@param lonCellIndex The longitudinal index of the cell to run stock ecology for
@param workingGridCellStocks A copy of the cohorts in the current grid cell
 * NB - acting on a copy here?
*/
void RunWithinCellStockEcology(unsigned latCellIndex, unsigned lonCellIndex, GridCellStockHandler workingGridCellStocks, int cellIndex)
        {
//            // Create a local instance of the stock ecology class
            EcologyStock MadingleyEcologyStock;
//            
//            // Initialise stock ecology - not needed!
//            MadingleyEcologyStock.InitializeEcology();
//
//            //The location of the acting stock
             vector<int> ActingStock(2);
//
//            // Get the list of functional group indices for autotroph stocks
            vector<int> AutotrophStockFunctionalGroups = StockFunctionalGroupDefinitions.GetFunctionalGroupIndex("Heterotroph/Autotroph", "Autotroph", false);
//
//            // Loop over autotroph functional groups
            for (int FunctionalGroup: AutotrophStockFunctionalGroups)
            {
                for (int ll = 0; ll < workingGridCellStocks[FunctionalGroup].size(); ll++)
                {
//                    // Get the position of the acting stock
                    ActingStock[0] = FunctionalGroup;
                    ActingStock[1] = ll;
//
//                    // Run stock ecology

                    MadingleyEcologyStock.RunWithinCellEcology(workingGridCellStocks, ActingStock, EcosystemModelGrid.GetCellEnvironment(
                        latCellIndex, lonCellIndex), EnvironmentalDataUnits, HumanNPPExtraction, StockFunctionalGroupDefinitions,
                        CurrentTimeStep, GlobalModelTimeStepUnit, ProcessTrackers[cellIndex].TrackProcesses, ProcessTrackers[cellIndex],TrackGlobalProcesses, CurrentMonth,
                        InitialisationFileStrings["OutputDetail"],SpecificLocations);
//

                    //workingGridCellStocks[ActingStock].TotalBiomass *= 0.75;
//
                }
            }

        }
//working with *copies* of the cohorts/stocks -why?
void RunWithinCellCohortEcology(unsigned latCellIndex, unsigned lonCellIndex, ThreadLockedParallelVariables& partial, 
            GridCellCohortHandler workingGridCellCohorts, GridCellStockHandler workingGridCellStocks,string outputDetail, int cellIndex, MadingleyModelInitialisation& initialisation)
        {
//
//
//            // Local instances of classes
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

            for (int i = 0; i < workingGridCellCohorts.size(); i++)
            {
                OriginalGridCellCohortsNumbers[i] = workingGridCellCohorts[i].size();
            }

            // Initialize ecology for stocks and cohorts
            MadingleyEcologyCohort.InitializeEcology(EcosystemModelGrid.GetCellEnvironment(latCellIndex, lonCellIndex)["Cell Area"][0],
                GlobalModelTimeStepUnit, DrawRandomly);

            // Create a jagged array indexed by functional groups to hold cohort indices
            CohortIndices.resize(CohortFunctionalGroupDefinitions.GetNumberOfFunctionalGroups());

//            // Loop over functional groups
            for (int ll = 0; ll < CohortFunctionalGroupDefinitions.GetNumberOfFunctionalGroups(); ll++)
            {
//                // Dimension the number of columns in each row of the jagged array to equal number of gridCellCohorts in each functional group
                if (workingGridCellCohorts[ll].size()==0)
                {
                    CohortIndices[ll].push_back(0); //correct??
                }
                else
                {
                    CohortIndices[ll].resize(workingGridCellCohorts[ll].size());
                }
                // Loop over gridCellCohorts in the functional group
                for (int kk = 0; kk < CohortIndices[ll].size(); kk++)
                {
                    // Fill jagged array with indices for each cohort
                    CohortIndices[ll][kk] = TotalCohortNumber;
                    TotalCohortNumber += 1;
                }

            }

            if (DrawRandomly)
            {
                // Randomly order the cohort indices
                RandomCohortOrder = Utilities.RandomlyOrderedCohorts(TotalCohortNumber);
            }
            else
            {
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
            for (int ll = 0; ll < RandomCohortOrder.size(); ll++)
            {

                // Locate the randomly chosen cohort within the array of lists of gridCellCohorts in the grid cell
                ActingCohort = Utilities.FindJaggedArrayIndex(RandomCohortOrder[ll], CohortIndices, TotalCohortNumber);

                // Perform all biological functions except dispersal (which is cross grid cell)
                if (workingGridCellCohorts[ActingCohort[0]].size()!=0 && workingGridCellCohorts[ActingCohort].CohortAbundance>ExtinctionThreshold)
                {
                    // Calculate number of cohorts in this functional group in this grid cell before running ecology
                    EcosystemModelParallelTempval1 = workingGridCellCohorts[ActingCohort[0]].size();

                    CohortActivity.AssignProportionTimeActive(workingGridCellCohorts[ActingCohort], EcosystemModelGrid.GetCellEnvironment(latCellIndex, lonCellIndex), CohortFunctionalGroupDefinitions, CurrentTimeStep, CurrentMonth);

//                    // Run ecology
                    MadingleyEcologyCohort.RunWithinCellEcology(workingGridCellCohorts, workingGridCellStocks,
                        ActingCohort, EcosystemModelGrid.GetCellEnvironment(latCellIndex, lonCellIndex),
                        EcosystemModelGrid.GetCellDeltas(latCellIndex, lonCellIndex),
                        CohortFunctionalGroupDefinitions, StockFunctionalGroupDefinitions, CurrentTimeStep,
                        ProcessTrackers[cellIndex],  partial, SpecificLocations,outputDetail, CurrentMonth, initialisation);

//                    // Update the properties of the acting cohort
                    MadingleyEcologyCohort.UpdateEcology(workingGridCellCohorts, workingGridCellStocks, ActingCohort,
                        EcosystemModelGrid.GetCellEnvironment(latCellIndex, lonCellIndex), EcosystemModelGrid.GetCellDeltas(
                        latCellIndex, lonCellIndex), CohortFunctionalGroupDefinitions, StockFunctionalGroupDefinitions, CurrentTimeStep,
                        ProcessTrackers[cellIndex]);

//                    // Add newly produced cohorts to the tracking variable
                    EcosystemModelParallelTempval2 += workingGridCellCohorts[ActingCohort[0]].size() - EcosystemModelParallelTempval1;
                    

//                    // Check that the mass of individuals in this cohort is still >= 0 after running ecology
                    assert(workingGridCellCohorts[ActingCohort].IndividualBodyMass() >= 0.0 && "Biomass < 0 for this cohort");
                }

                // Check that the mass of individuals in this cohort is still >= 0 after running ecology
                if (workingGridCellCohorts[ActingCohort[0]].size()>0 )assert(workingGridCellCohorts[ActingCohort].IndividualBodyMass() >= 0.0 && "Biomass < 0 for this cohort");
            }


            // Update diagnostics of productions
            partial.Productions += EcosystemModelParallelTempval2;

            RunExtinction(latCellIndex, lonCellIndex, partial, workingGridCellCohorts, cellIndex);



            // Merge cohorts, if necessary
            if (workingGridCellCohorts.GetNumberOfCohorts() > initialisation.MaxNumberOfCohorts)
            {
                partial.Combinations = CohortMerger.MergeToReachThresholdFast(workingGridCellCohorts, workingGridCellCohorts.GetNumberOfCohorts(), initialisation.MaxNumberOfCohorts);

                //Run extinction a second time to remove those cohorts that have been set to zero abundance when merging
                RunExtinction(latCellIndex, lonCellIndex, partial, workingGridCellCohorts, cellIndex);
            }
            else
                partial.Combinations = 0;
            
            // Write out the updated cohort numbers after all ecological processes have occured
            EcosystemModelGrid.SetGridCellCohorts(workingGridCellCohorts, latCellIndex, lonCellIndex);
        }

/** \brief Carries out extinction on cohorts that have an abundance below a defined extinction threshold */
void RunExtinction(unsigned latCellIndex, unsigned lonCellIndex, ThreadLockedParallelVariables& partial,
            GridCellCohortHandler workingGridCellCohorts, int cellIndex)
        {
            bool VarExists;

            // Loop over cohorts and remove any whose abundance is below the extinction threshold
            for (int kk = 0; kk < CohortFunctionalGroupDefinitions.GetNumberOfFunctionalGroups(); kk++)
            {
                // Create a list to hold the cohorts to remove
                vector<int> CohortIndicesToRemove;

                // Loop through each cohort in the functional group
                for (int ll = 0; ll < workingGridCellCohorts[kk].size(); ll++)
                {
                    // If cohort abundance is less than the extinction threshold then add to the list for extinction
                    if (workingGridCellCohorts[kk][ll].CohortAbundance> ExtinctionThreshold || workingGridCellCohorts[kk][ll].IndividualBodyMass() <= 1.e-300)
                    {
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
                for (int ll = (CohortIndicesToRemove.size() - 1); ll >= 0; ll--)
                {
//                    // Add biomass of the extinct cohort to the organic matter pool
                    EcosystemModelGrid.SetEnviroLayer("Organic Pool", 0, EcosystemModelGrid.GetEnviroLayer("Organic Pool", 0, latCellIndex, lonCellIndex, VarExists) +
                        (workingGridCellCohorts[kk][CohortIndicesToRemove[ll]].IndividualBodyMass() + workingGridCellCohorts[kk][CohortIndicesToRemove[ll]].IndividualReproductivePotentialMass) * workingGridCellCohorts[kk][CohortIndicesToRemove[ll]].CohortAbundance, latCellIndex, lonCellIndex);
                    assert(EcosystemModelGrid.GetEnviroLayer("Organic Pool", 0, latCellIndex, lonCellIndex, VarExists) >= 0 && "Organic pool < 0");
//
//                    if (ProcessTrackers[cellIndex].TrackProcesses && SpecificLocations == true)
//                        ProcessTrackers[cellIndex].RecordExtinction(latCellIndex, lonCellIndex, CurrentTimeStep, workingGridCellCohorts[kk][CohortIndicesToRemove[ll]].Merged, workingGridCellCohorts[kk][CohortIndicesToRemove[ll]].CohortID);

                    // Remove the extinct cohort from the list of cohorts
                    workingGridCellCohorts[kk].erase(workingGridCellCohorts[kk].begin()+CohortIndicesToRemove[ll]);


                }

            }

        }

void RunWithinCellDispersalOnly(unsigned latCellIndex, unsigned lonCellIndex, ThreadLockedParallelVariables& partial,
                   GridCellCohortHandler workingGridCellCohorts, GridCellStockHandler workingGridCellStocks)
        {
//            // Merge cohorts. Requires cohorts to be identical, for testing purposes (remember that they don't grow etc)
//            // SHOULD ONLY BE RUN FOR RESPONSIVE DISPERSAL TESTING
//            //partial.Combinations = Merger.MergeForResponsiveDispersalOnly(workingGridCellCohorts);
//
//            // Loop over cohorts and remove any whose abundance is below the extinction threshold
//            for (int kk = 0; kk < CohortFunctionalGroupDefinitions.GetNumberOfFunctionalGroups(); kk++)
//            {
//                // Create a list to hold the cohorts to remove
//                List<int> CohortIndicesToRemove = new List<int>();
//
//                // Loop through each cohort in the functional group
//                for (int ll = 0; ll < workingGridCellCohorts[kk].Count; ll++)
//                {
//                    // If cohort abundance is less than the extinction threshold then add to the list for extinction
//                    if (workingGridCellCohorts[kk][ll].CohortAbundance <= _ExtinctionThreshold)
//                    {
//                        CohortIndicesToRemove.Add(ll);
//
//                        partial.Extinctions += 1;
//                    }
//                }
//
//                // Note that we don't keep track of the organic biomass pool if running dispersal only, since there are cohorts with strange biomasses
//                for (int ll = (CohortIndicesToRemove.Count - 1); ll >= 0; ll--)
//                {
//                    // Remove the extinct cohort from the list of cohorts
//                    workingGridCellCohorts[kk].RemoveAt(CohortIndicesToRemove[ll]);
//                }
//
//            }
//
//
//            // Write out the updated cohort numbers after all ecological processes have occured
//            EcosystemModelGrid.SetGridCellCohorts(workingGridCellCohorts, latCellIndex, lonCellIndex);
        }

/** \brief Run ecological processes that operate across grid cells */
void RunCrossGridCellEcology(unsigned& dispersals, bool dispersalOnly, MadingleyModelInitialisation& modelInitialisation)
        {
//            // If we are running specific locations, then we do not run dispersal
//            if (SpecificLocations != true)
//            {
//                if (RunGridCellsInParallel)
//                {
//                    // Loop through each grid cell, and run dispersal for each.
//                    // Note that currently dispersal is not parallelised, although it could be (though care would need to be taken to ensure that necessary variables are thread-locked
//                    //for (int ii = 0; ii < _CellList.Count; ii++)
//                    Parallel.For(0, _CellList.Count, ii =>
//                    {
//                        EcologyCrossGridCell TempMadingleyEcologyCrossGridCell = new EcologyCrossGridCell();
//                        // Initialise cross grid cell ecology
//                        TempMadingleyEcologyCrossGridCell.InitializeCrossGridCellEcology(_GlobalModelTimeStepUnit, DrawRandomly, modelInitialisation);
//
//                        //Initialise the delta for dispersal lists for this grid cell
//                        EcosystemModelGrid.DeltaFunctionalGroupDispersalArray[_CellList[ii][0], _CellList[ii][1]] = new List<uint>();
//                        EcosystemModelGrid.DeltaCohortNumberDispersalArray[_CellList[ii][0], _CellList[ii][1]] = new List<uint>();
//                        EcosystemModelGrid.DeltaCellToDisperseToArray[_CellList[ii][0], _CellList[ii][1]] = new List<uint[]>();
//
//                        // We have looped through individal cells and calculated ecological processes for each. Now do this for cross grid cell processes
//                        TempMadingleyEcologyCrossGridCell.RunCrossGridCellEcology(_CellList[ii], EcosystemModelGrid, dispersalOnly,
//                            CohortFunctionalGroupDefinitions, StockFunctionalGroupDefinitions, CurrentMonth);
//                    });
//                }
//                else
//                {
//                    // Loop through each grid cell, and run dispersal for each.
//                    // Note that currently dispersal is not parallelised, although it could be (though care would need to be taken to ensure that necessary variables are thread-locked
    for (int ii = 0; ii < CellList.size(); ii++)
                    {
                        //Initialise the delta for dispersal lists for this grid cell
                        
                       // EcosystemModelGrid.DeltaFunctionalGroupDispersalArray[CellList[ii][0], CellList[ii][1]].empty();
                       // EcosystemModelGrid.DeltaCohortNumberDispersalArray[CellList[ii][0], CellList[ii][1]].empty() ;
                       // EcosystemModelGrid.DeltaCellToDisperseToArray[CellList[ii][0], CellList[ii][1]].empty() ;
                        

                        // We have looped through individal cells and calculated ecological processes for each. Now do this for cross grid cell processes
                        MadingleyEcologyCrossGridCell.RunCrossGridCellEcology(CellList[ii], EcosystemModelGrid, dispersalOnly,
                            CohortFunctionalGroupDefinitions, StockFunctionalGroupDefinitions, CurrentMonth);
                    }
//                }
//                // Apply the changes in the delta arrays from dispersal
                MadingleyEcologyCrossGridCell.UpdateCrossGridCellEcology(EcosystemModelGrid, dispersals);
                            

                //            }
//
        }

/** \brief Make a record of the properties of the intial model cohorts in the new cohorts output file */
void RecordInitialCohorts()
        {
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

};
#endif