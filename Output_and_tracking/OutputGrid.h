#ifndef OUTPUTGRID_H
#define OUTPUTGRID_H
/** \file OutputGrid.h
 * \brief the OutputGrid header file
 */






//

//


//

//
//
//namespace Madingley
//{
//
//
//    class OutputGrid
//    {
/** \brief
Designates the level of output detail
*/
//        private enum OutputDetailLevel { Low, Medium, High };
//
/** \brief
An instance of the enumerator to designate output detail level
*/
//        OutputDetailLevel ModelOutputDetail;
//
/** \brief
A dataset to store the live screen view
*/
//        private DataSet DataSetToViewLive;
//
/** \brief
The path to the output folder
*/
//        private string _OutputPath;
/** \brief
Get the path to the output folder
*/
//         string OutputPath { get { return _OutputPath; } }
//
/** \brief
Dataset object to handle grid-based outputs
*/
//        private DataSet GridOutput;
//
/** \brief
The cohort traits to be considered in the outputs
*/
//        private string[] CohortTraits;
//
/** \brief
All unique values of the traits to be considered in outputs
*/
//        private SortedDictionary<string, string[]> CohortTraitValues;
//
/** \brief
The stock traits to be considered in the outputs
*/
//        private string[] StockTraits;
//
/** \brief
All unique values of the traits to be considered in the outputs
*/
//        private SortedDictionary<string, string[]> StockTraitValues;
//
/** \brief
Holds a list of the functional group indices correpsonding to each unique cohort trait
*/
//        private map<string, int[]> CohortTraitIndices = new SortedList<string, int[]>();
//
/** \brief
Holds a list of the functional group indices corresponding to each unique stock trait
*/
//        private map<string, int[]> StockTraitIndices = new SortedList<string, int[]>();
//
/** \brief
Grid of (log) total biomass densities of cohorts in individual grid cells
*/
//        private double[,] LogBiomassDensityGridCohorts;
//
/** \brief
Grid of (log) total biomass densities of stocks in individual grid cells
*/
//        private double[,] LogBiomassDensityGridStocks;
//
/** \brief
Grid of (log) total biomass densities of both stocks and cohorts in individual grid cells
*/
//        private double[,] LogBiomassDensityGrid;
//
/** \brief
Grod of (log) abundance densities of cohorts in individual grid cells
*/
//        private double[,] LogAbundanceDensityGridCohorts;
//
/** \brief
Grids of total biomass densities in individual grid cells, arranged by trait value
*/
//        private map<string, double[,]> BiomassDensityGrid = new SortedList<string, double[,]>();
//
/** \brief
Grids of total densities in individual grid cells, arranged by trait value
*/
//        private map<string, double[,]> AbundanceDensityGrid = new SortedList<string, double[,]>();
//
//        private double[,] Realm;
//
//        private double[,] FrostDays;
//
//        private double[,] FracEvergreen;
//        
/** \brief
The time steps in this model simulation
*/
//        private float[] TimeSteps;
//
/** \brief
An instance of the class to convert data between arrays and SDS objects
*/
//        private ArraySDSConvert DataConverter;
//
/** \brief
Instance of the class to create SDS objects
*/
//        private CreateSDSObject SDSCreator;
//
/** \brief
An instance of the class for viewing grid results
*/
//        private ViewGrid GridViewer;
//
/** \brief
Whether to display live outputs during this model run
*/
//        private Boolean LiveOutputs;
//
//         OutputGrid(string outputDetail, MadingleyModelInitialisation modelInitialisation)
//        {
//            // Set the output path
//            _OutputPath = modelInitialisation.OutputPath;
//
//            // Initialise the grid viewer
//            GridViewer = new ViewGrid();
//
//            // Set the output detail level
//            if (outputDetail == "low")
//                ModelOutputDetail = OutputDetailLevel.Low;
//            else if (outputDetail == "medium")
//                ModelOutputDetail = OutputDetailLevel.Medium;
//            else if (outputDetail == "high")
//                ModelOutputDetail = OutputDetailLevel.High;
//            else
//                Debug.Fail("Specified output detail level is not valid, must be 'low', 'medium' or 'high'");
//
//            // Initialise the data converter
//            DataConverter = new ArraySDSConvert();
//
//            // Initialise the SDS object creator
//            SDSCreator = new CreateSDSObject();
//
//            // Set the local variable designating whether to display live outputs during this model run
//            if (modelInitialisation.LiveOutputs)
//                LiveOutputs = true;
//
//        }
//
/** \brief
Spawn dataset viewer for the live outputs
*/
//         void SpawnDatasetViewer()
//        {
//            Console.WriteLine("Spawning Dataset Viewer\n");
//
//            // Intialise the SDS object for the live view
//            DataSetToViewLive = SDSCreator.CreateSDSInMemory(true);
//
//            DataSetToViewLive.Metadata["VisualHints"] =
//                    "\"Biomass density\"(Longitude,Latitude) Style:Colormap; Palette:#000040=0,Blue=0.1,Green=0.2661,#FA8000=0.76395,Red; Transparency:0.38";
//
//            // Start viewing
//            GridViewer.AsynchronousView(ref DataSetToViewLive, "");
//
//        }
//
//         void SetupOutputs(ModelGrid ecosystemModelGrid, string outputFilesSuffix, unsigned numTimeSteps,
//            FunctionalGroupDefinitions cohortFunctionalGroupDefinitions, FunctionalGroupDefinitions 
//            stockFunctionalGroupDefinitions)
//        {
//            // Create an SDS object to hold grid outputs
//            GridOutput = SDSCreator.CreateSDS("netCDF", "GridOutputs" + outputFilesSuffix, _OutputPath);
//
//            // Initilalise trait-based outputs
//            InitialiseTraitBasedOutputs(cohortFunctionalGroupDefinitions, stockFunctionalGroupDefinitions);
//
//            // Create vector to hold the values of the time dimension
//            TimeSteps = new float[numTimeSteps + 1];
//
//            // Set the first value to be 0 (this will hold initial outputs)
//            TimeSteps[0] = 0;
//
//            // Fill other values from 0 (this will hold outputs during the model run)
//            for (int i = 1; i < numTimeSteps + 1; i++)
//            {
//                TimeSteps[i] = i;
//            }
//
//            // Declare vectors for geographical dimension data
//            float[] outLats = new float[ecosystemModelGrid.NumLatCells];
//            float[] outLons = new float[ecosystemModelGrid.NumLonCells];
//
//            // Populate the dimension variable vectors with cell centre latitude and longitudes
//            for (int i = 0; i < ecosystemModelGrid.NumLatCells; i++)
//            {
//                outLats[i] = ecosystemModelGrid.Lats[i] + (ecosystemModelGrid.LatCellSize / 2);
//            }
//
//            for (int jj = 0; jj < ecosystemModelGrid.NumLonCells; jj++)
//            {
//                outLons[jj] = ecosystemModelGrid.Lons[jj] + (ecosystemModelGrid.LonCellSize / 2);
//            }
//
//            //GridOutputArray = new double[ecosystemModelGrid.NumLatCells, ecosystemModelGrid.NumLonCells, numTimeSteps];
//
//            // Add output variables that are dimensioned geographically and temporally to grid output file
//            string[] GeographicalDimensions = { "Latitude", "Longitude", "Time step" };
//            DataConverter.AddVariable(GridOutput, "Biomass density", 3, GeographicalDimensions, 0, outLats, outLons, TimeSteps);
//            DataConverter.AddVariable(GridOutput, "Abundance density", 3, GeographicalDimensions, 0, outLats, outLons, TimeSteps);
//
//            // Initialise the arrays that will be used for the grid-based outputs
//            LogBiomassDensityGridCohorts = new double[ecosystemModelGrid.NumLatCells, ecosystemModelGrid.NumLonCells];
//            LogBiomassDensityGridStocks = new double[ecosystemModelGrid.NumLatCells, ecosystemModelGrid.NumLonCells];
//            LogBiomassDensityGrid = new double[ecosystemModelGrid.NumLatCells, ecosystemModelGrid.NumLonCells];
//            LogAbundanceDensityGridCohorts = new double[ecosystemModelGrid.NumLatCells, ecosystemModelGrid.NumLonCells];
//            FracEvergreen = new double[ecosystemModelGrid.NumLatCells, ecosystemModelGrid.NumLonCells];
//            Realm = new double[ecosystemModelGrid.NumLatCells, ecosystemModelGrid.NumLonCells];
//
//            // Temporary outputs for checking plant model
//            DataConverter.AddVariable(GridOutput, "Fraction year frost", 3, GeographicalDimensions, ecosystemModelGrid.GlobalMissingValue, outLats, outLons, TimeSteps);
//            DataConverter.AddVariable(GridOutput, "Fraction evergreen", 3, GeographicalDimensions, ecosystemModelGrid.GlobalMissingValue, outLats, outLons, TimeSteps);
//            DataConverter.AddVariable(GridOutput, "Realm", 3, GeographicalDimensions, ecosystemModelGrid.GlobalMissingValue, outLats, outLons, TimeSteps);
//
//            // Set up outputs for medium or high output levels
//            if ((ModelOutputDetail == OutputDetailLevel.Medium) || (ModelOutputDetail == OutputDetailLevel.High))
//            {
//                double[,] temp = new double[ecosystemModelGrid.NumLatCells, ecosystemModelGrid.NumLonCells];
//                foreach (string TraitValue in CohortTraitIndices.Keys)
//                {
//                    BiomassDensityGrid.Add(TraitValue, temp);
//                    DataConverter.AddVariable(GridOutput, TraitValue + "biomass density", 3, GeographicalDimensions, 0, outLats, outLons, TimeSteps);
//                    AbundanceDensityGrid.Add(TraitValue, temp);
//                    DataConverter.AddVariable(GridOutput, TraitValue + "abundance density", 3, GeographicalDimensions, 0, outLats, outLons, TimeSteps);
//                }
//                foreach (string TraitValue in StockTraitIndices.Keys)
//                {
//                    BiomassDensityGrid.Add(TraitValue, temp);
//                    DataConverter.AddVariable(GridOutput, TraitValue + "biomass density", 3, GeographicalDimensions, 0, outLats, outLons, TimeSteps);
//                }
//            }
//        }
//
/** \brief
Set up the necessary architecture for generating outputs arranged by trait value
*/
@param cohortFunctionalGroupDefinitions Functional group definitions for cohorts in the model 
@param stockFunctionalGroupDefinitions Functional group definitions for stocks in the model 
//        private void InitialiseTraitBasedOutputs(FunctionalGroupDefinitions cohortFunctionalGroupDefinitions, FunctionalGroupDefinitions
//            stockFunctionalGroupDefinitions)
//        {
//            // Define the cohort traits that will be used to separate outputs
//            CohortTraits = new string[3] { "Nutrition source", "Endo/Ectotherm", "Reproductive strategy" };
//
//            // Declare a sorted dictionary to hold all unique trait values
//            CohortTraitValues = new SortedDictionary<string, string[]>();
//
//            // Add all unique trait values to the sorted dictionary
//            foreach (string Trait in CohortTraits)
//            {
//                CohortTraitValues.Add(Trait, cohortFunctionalGroupDefinitions.GetUniqueTraitValues(Trait));
//
//            }
//
//            // Get the list of functional group indices corresponding to each unique trait value
//            foreach (string Trait in CohortTraits)
//            {
//                foreach (string TraitValue in CohortTraitValues[Trait])
//                {
//                    CohortTraitIndices.Add(TraitValue, cohortFunctionalGroupDefinitions.GetFunctionalGroupIndex(Trait, TraitValue, false));
//                }
//            }
//
//            // Define the stock traits that will be used to separate outputs
//            StockTraits = new string[2] { "Heterotroph/Autotroph" ,"Leaf strategy"};
//
//            // Re-initialise the sorted dictionary to hold all unique trait values
//            StockTraitValues = new SortedDictionary<string, string[]>();
//
//            // Add all unique stock trait values to the sorted dictionary
//            foreach (string Trait in StockTraits)
//            {
//                StockTraitValues.Add(Trait, stockFunctionalGroupDefinitions.GetUniqueTraitValues(Trait));
//            }
//
//            // Get the list of functional group indices corresponding to each unique trait value
//            foreach (string Trait in StockTraits)
//            {
//                foreach (string TraitValue in StockTraitValues[Trait])
//                {
//                    StockTraitIndices.Add(TraitValue, stockFunctionalGroupDefinitions.GetFunctionalGroupIndex(Trait, TraitValue, false));
//                }
//            }
//
//        }
//
//        private void CalculateOutputs(ModelGrid ecosystemModelGrid, FunctionalGroupDefinitions cohortFunctionalGroupDefinitions,
//            FunctionalGroupDefinitions stockFunctionalGroupDefinitions, List<unsigned[]> cellIndices, MadingleyModelInitialisation initialisation)
//        {
//            // Get grids of the total biomass densities of all stocks and all cohorts in each grid cell
//            LogBiomassDensityGridCohorts = ecosystemModelGrid.GetStateVariableGridLogDensityPerSqKm("Biomass", "NA", cohortFunctionalGroupDefinitions.
//                AllFunctionalGroupsIndex, cellIndices, "cohort", initialisation);
//            LogBiomassDensityGridStocks = ecosystemModelGrid.GetStateVariableGridLogDensityPerSqKm("Biomass", "NA", stockFunctionalGroupDefinitions.
//                AllFunctionalGroupsIndex, cellIndices, "stock", initialisation);
//
//            // Get grids of total abundance densities of all stocks and all cohorts in each grid cell
//            LogAbundanceDensityGridCohorts = ecosystemModelGrid.GetStateVariableGridLogDensityPerSqKm("Abundance", "NA", cohortFunctionalGroupDefinitions.
//                 AllFunctionalGroupsIndex, cellIndices, "cohort", initialisation);
//
//            // Loop over grid cells and add stock and cohort biomass density to get the total of all biomass densities
//            for (int ii = 0; ii < ecosystemModelGrid.NumLatCells; ii++)
//            {
//                for (int jj = 0; jj < ecosystemModelGrid.NumLonCells; jj++)
//                {
//                    LogBiomassDensityGrid[ii, jj] = Math.Log(Math.Exp(LogBiomassDensityGridCohorts[ii, jj]) + Math.Exp(LogBiomassDensityGridStocks[ii, jj]));
//                }
//            }
//
//            string[] Keys = CohortTraitIndices.Keys.ToArray();
//            foreach (string Key in Keys)
//            {
//                BiomassDensityGrid[Key] = ecosystemModelGrid.GetStateVariableGridLogDensityPerSqKm("Biomass", "NA", CohortTraitIndices[Key], cellIndices, "cohort", initialisation);
//                AbundanceDensityGrid[Key] = ecosystemModelGrid.GetStateVariableGridLogDensityPerSqKm("Abundance", "NA", CohortTraitIndices[Key], cellIndices, "cohort", initialisation);
//            }
//            
//            Keys = StockTraitIndices.Keys.ToArray();
//            foreach (string Key in Keys)
//            {
//                BiomassDensityGrid[Key] = ecosystemModelGrid.GetStateVariableGridLogDensityPerSqKm("Biomass", "NA", StockTraitIndices[Key], cellIndices, "stock", initialisation);
//            }
//
//            
//            // Temporary outputs to check plant model
//
//            Realm = ecosystemModelGrid.GetEnviroGrid("Realm", 0);
//
//            FrostDays = ecosystemModelGrid.GetEnviroGrid("Fraction Year Frost", 0);
//
//            for (int i = 0; i < ecosystemModelGrid.NumLatCells; i++)
//            {
//                for (int j = 0; j < ecosystemModelGrid.NumLonCells; j++)
//                {
//                    FracEvergreen[i, j] = BiomassDensityGrid["evergreen"][i, j] / BiomassDensityGrid["autotroph"][i, j];
//                }
//            }
//
//            Realm = ecosystemModelGrid.GetEnviroGrid("Realm", 0);
//
//        }
//
//         void InitialOutputs(ModelGrid ecosystemModelGrid, FunctionalGroupDefinitions cohortFunctionalGroupDefinitions, FunctionalGroupDefinitions
//            stockFunctionalGroupDefinitions, List<unsigned[]> cellIndices, MadingleyModelInitialisation initialisation)
//        {
//            Console.WriteLine("Writing initial grid outputs...");
//
//            // Calculate the output variables
//            CalculateOutputs(ecosystemModelGrid, cohortFunctionalGroupDefinitions, stockFunctionalGroupDefinitions, cellIndices, initialisation);
//
//            // Write the total biomass of cohorts to the live display
//            if (LiveOutputs)
//            {
//                DataConverter.Array2DToSDS2D(LogBiomassDensityGridCohorts, "Log(Biomass density, g/km^2)", ecosystemModelGrid.Lats,
//                    ecosystemModelGrid.Lons, ecosystemModelGrid.GlobalMissingValue, DataSetToViewLive);
//            }
//
//            // Add the grid of total biomass in cells to the file dataset
//            DataConverter.Array2DToSDS3D(LogBiomassDensityGridCohorts, "Biomass density", new string[] { "Latitude", "Longitude", "Time step" }, 0, 
//                ecosystemModelGrid.GlobalMissingValue, GridOutput);
//
//            DataConverter.Array2DToSDS3D(LogAbundanceDensityGridCohorts, "Abundance density", new string[] { "Latitude", "Longitude", "Time step" }, 0,
//                ecosystemModelGrid.GlobalMissingValue, GridOutput);
//
//            // Temporary outputs to check plant model
//            DataConverter.Array2DToSDS3D(Realm, "Realm", new string[] { "Latitude", "Longitude", "Time step" },
//                0, ecosystemModelGrid.GlobalMissingValue, GridOutput);
//
//            // File outputs for medium and high detail levels
//            if ((ModelOutputDetail == OutputDetailLevel.Medium) || (ModelOutputDetail == OutputDetailLevel.High))
//            {
//                // Add the biomass grids for individual trait combinations to the file dataset
//                foreach (var Key in BiomassDensityGrid.Keys)
//                {
//                    DataConverter.Array2DToSDS3D(BiomassDensityGrid[Key], Key + "biomass density", new string[]
//                    { "Latitude", "Longitude", "Time step" }, 0, ecosystemModelGrid.GlobalMissingValue, GridOutput);
//                }
//
//                // Add the abundance density grid
//                foreach (var Key in AbundanceDensityGrid.Keys)
//                {
//
//                    DataConverter.Array2DToSDS3D(AbundanceDensityGrid[Key], Key + "abundance density",
//                        new string[] { "Latitude", "Longitude", "Time step" },
//                        0,
//                        ecosystemModelGrid.GlobalMissingValue,
//                        GridOutput);
//                }
//            }
//           
//        }
//
//        
//         void TimeStepOutputs(ModelGrid ecosystemModelGrid, FunctionalGroupDefinitions cohortFunctionalGroupDefinitions, FunctionalGroupDefinitions
//            stockFunctionalGroupDefinitions, List<unsigned[]> cellIndices, unsigned currentTimeStep, MadingleyModelInitialisation initialisation)
//        {
//            // Calculate the output variables for this time step
//            CalculateOutputs(ecosystemModelGrid, cohortFunctionalGroupDefinitions, stockFunctionalGroupDefinitions, cellIndices, initialisation);
//
//            // Write the total biomass of cohorts to the live display
//            if (LiveOutputs)
//            {
//                DataConverter.Array2DToSDS2D(LogBiomassDensityGridCohorts, "Log(Biomass density, g/km^2)", ecosystemModelGrid.Lats,
//                    ecosystemModelGrid.Lons, 0, DataSetToViewLive);
//            }
//
//            Console.WriteLine("Writing grid ouputs to file...\n");
//
//            // Add the grid of total biomass in cells to the file dataset
//            DataConverter.Array2DToSDS3D(LogBiomassDensityGridCohorts, "Biomass density", new string[] { "Latitude", "Longitude", "Time step" },
//                (int)currentTimeStep + 1, 0, GridOutput);
//
//            // Add the grid of total abudance in cells to the file dataset
//            DataConverter.Array2DToSDS3D(LogAbundanceDensityGridCohorts, "Abundance density", new string[] { "Latitude", "Longitude", "Time step" },
//                (int)currentTimeStep + 1, 0, GridOutput);
//
//            // Temporary outputs to check plant model
//            DataConverter.Array2DToSDS3D(Realm, "Realm", new string[] { "Latitude", "Longitude", "Time step" },
//                (int)currentTimeStep + 1, ecosystemModelGrid.GlobalMissingValue, GridOutput);
//            DataConverter.Array2DToSDS3D(FrostDays, "Fraction year frost", new string[] { "Latitude", "Longitude", "Time step" },
//                (int)currentTimeStep + 1, ecosystemModelGrid.GlobalMissingValue, GridOutput);
//            DataConverter.Array2DToSDS3D(FracEvergreen, "Fraction evergreen", new string[] { "Latitude", "Longitude", "Time step" },
//                (int)currentTimeStep + 1, ecosystemModelGrid.GlobalMissingValue, GridOutput);
//            DataConverter.Array2DToSDS3D(Realm, "Realm", new string[] { "Latitude", "Longitude", "Time step" },
//                (int)currentTimeStep + 1, ecosystemModelGrid.GlobalMissingValue, GridOutput);
//
//
//
//            if ((ModelOutputDetail == OutputDetailLevel.Medium) || (ModelOutputDetail == OutputDetailLevel.High))
//            {
//                foreach (var TraitValue in BiomassDensityGrid.Keys)
//                {
//                    // Add the biomass grids for individual trait combinations to the file dataset
//                    DataConverter.Array2DToSDS3D(BiomassDensityGrid[TraitValue], TraitValue + "biomass density", new string[] 
//                    { "Latitude", "Longitude", "Time step" }, (int)currentTimeStep+1, ecosystemModelGrid.GlobalMissingValue, GridOutput);
//                }
//
//                foreach (var Key in AbundanceDensityGrid.Keys)
//                {
//
//                    DataConverter.Array2DToSDS3D(AbundanceDensityGrid[Key], Key + "abundance density",
//                        new string[] { "Latitude", "Longitude", "Time step" },
//                        (int)currentTimeStep + 1,
//                        ecosystemModelGrid.GlobalMissingValue,
//                        GridOutput);
//                }
//            }
//
//
//        }
//
//         void FinalOutputs()
//        {
//            // Dispose of the grid outputs dataset
//            GridOutput.Dispose();
//        }
//
//    }
//}
#endif
