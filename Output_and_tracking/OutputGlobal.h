#ifndef OUTPUTGLOBAL_H
#define OUTPUTGLOBAL_H
/** \file OutputGlobal.h
 * \brief the OutputGlobal header file
 */





//


//


//
//namespace Madingley
//{
//    class OutputGlobal
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
A dataset to store the basic outputs to file
*/
//        private DataSet BasicOutput;
//        
/** \brief
The time steps in this model simulation
*/
//        private float[] TimeSteps; 
//
/** \brief
The path to the output folder
*/
//        private string _OutputPath;
/** \brief
Get the path to the output folder
*/
//         string OutputPath
//        { get { return _OutputPath; } }
//
/** \brief
The suffix to apply to all outputs from this grid cell
*/
//        private string _OutputSuffix;
/** \brief
Get the suffix for ouputs for this grid cell
*/
//         string OutputSuffix
//        { get { return _OutputSuffix; } }
//
/** \brief
The total abundance across all cohorts in the grid cell
*/
//        private double TotalAbundance;
//
/** \brief
The total living biomass of all cohorts and stock in the grid cell
*/
//        private double TotalLivingBiomass;
//
/** \brief
The total biomass, living and non-living in the grid cell
*/
//        private double TotalBiomass;
//
/** \brief
The total biomass in the organic pool
*/
//        private double OrganicPoolOut;
//
/** \brief
The total biomass in the respiratory pool
*/
//        private double RespiratoryPoolOut;
//
/** \brief
Total number of cohorts in the grid cell
*/
//        private double TotalNumberOfCohorts;
//
/** \brief
Total number of stocks in the grid cell
*/
//        private double TotalNumberOfStocks;
//
/** \brief
Number of cohorts that became extinct in a time step
*/
//        private double NumberOfCohortsExtinct;
//
/** \brief
Number of cohorts produced in a time step
*/
//        private double NumberOfCohortsProduced;
//
/** \brief
Number of cohorts merged in a time step
*/
//        private double NumberOfCohortsCombined;
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
Constructor for the global output class
*/
@param outputDetail The level of detail to be used in model outputs 
@param modelInitialisation Model intialisation object 
//         OutputGlobal(string outputDetail, MadingleyModelInitialisation modelInitialisation)
//        {
//            // Set the output path
//            _OutputPath = modelInitialisation.OutputPath;
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
//        }
//
/** \brief
Set up all outputs (live, console and file) prior to the model run
*/
@param numTimeSteps The number of time steps in the model run 
@param ecosystemModelGrid The model grid 
@param outputFilesSuffix The suffix to be applied to all output files 
//         void SetupOutputs(unsigned numTimeSteps,ModelGrid ecosystemModelGrid, string outputFilesSuffix)
//        {
//            Console.WriteLine("Setting up global outputs...\n");
//
//            // Set the suffix for all output files
//            _OutputSuffix = outputFilesSuffix + "_Global";
//
//            // Create an SDS object to hold total abundance and biomass data
//            BasicOutput = SDSCreator.CreateSDS("netCDF", "BasicOutputs" + _OutputSuffix, _OutputPath);
//
//            // Create vector to hold the values of the time dimension
//            TimeSteps = new float[numTimeSteps + 1];
//
//            // Set the first value to be -1 (this will hold initial outputs)
//            TimeSteps[0] = 0;
//
//            // Fill other values from 0 (this will hold outputs during the model run)
//            for (int i = 1; i < numTimeSteps + 1; i++)
//            {
//                TimeSteps[i] = i;
//            }
//
//            // Add basic variables to the output file, dimensioned by time
//            string[] TimeDimension = { "Time step" };
//            DataConverter.AddVariable(BasicOutput, "Total living biomass", "Kg / km^2", 1, TimeDimension, ecosystemModelGrid.GlobalMissingValue, TimeSteps);
//            DataConverter.AddVariable(BasicOutput, "Organic matter pool", "Kg / km^2", 1, TimeDimension, ecosystemModelGrid.GlobalMissingValue, TimeSteps);
//            DataConverter.AddVariable(BasicOutput, "Respiratory CO2 pool", "Kg / km^2", 1, TimeDimension, ecosystemModelGrid.GlobalMissingValue, TimeSteps);
//            DataConverter.AddVariable(BasicOutput, "Number of cohorts extinct", "", 1, TimeDimension, ecosystemModelGrid.GlobalMissingValue, TimeSteps);
//            DataConverter.AddVariable(BasicOutput, "Number of cohorts produced", "", 1, TimeDimension, ecosystemModelGrid.GlobalMissingValue, TimeSteps);
//            DataConverter.AddVariable(BasicOutput, "Number of cohorts combined", "", 1, TimeDimension, ecosystemModelGrid.GlobalMissingValue, TimeSteps);
//            DataConverter.AddVariable(BasicOutput, "Number of cohorts in model", "", 1, TimeDimension, ecosystemModelGrid.GlobalMissingValue, TimeSteps);
//            DataConverter.AddVariable(BasicOutput, "Number of stocks in model", "", 1, TimeDimension, ecosystemModelGrid.GlobalMissingValue, TimeSteps);
//
//        }
//
/** \brief
Calculates the variables to output
*/
@param cohortFunctionalGroupDefinitions The functional group definitions of cohorts in the model 
@param stockFunctionalGroupDefinitions The functional group definitions of stocks in the model 
@param ecosystemModelGrid The model grid 
@param cellIndices The list of indices of active cells in the model grid 
@param globalDiagnosticVariables Global diagnostic variables 
//         void CalculateOutputs(FunctionalGroupDefinitions cohortFunctionalGroupDefinitions, FunctionalGroupDefinitions
//            stockFunctionalGroupDefinitions, ModelGrid ecosystemModelGrid, List<unsigned[]> cellIndices, map<string,double>
//            globalDiagnosticVariables, MadingleyModelInitialisation initialisation)
//        {
//            // Get all cohort functional group indices in the model
//            int[] CohortFunctionalGroupIndices = cohortFunctionalGroupDefinitions.AllFunctionalGroupsIndex;
//
//            // Get all stock functional group indices in the model
//            int[] StockFunctionalGroupIndices = stockFunctionalGroupDefinitions.AllFunctionalGroupsIndex;
//
//            // Reset total abundance, biomass and pool biomasses
//            TotalAbundance = 0.0;
//            TotalBiomass = 0.0;
//            TotalLivingBiomass = 0.0;
//            OrganicPoolOut = 0.0;
//            RespiratoryPoolOut = 0.0;
//
//            // Add total cohort biomass and total stock biomass to the total biomass tracker
//            TotalLivingBiomass += ecosystemModelGrid.StateVariableGridTotal("Biomass", "NA", CohortFunctionalGroupIndices, cellIndices, "cohort", initialisation);
//            TotalLivingBiomass += ecosystemModelGrid.StateVariableGridTotal("Biomass", "NA", StockFunctionalGroupIndices, cellIndices, "stock", initialisation);
//
//            // Add total cohort abundance to the total abundance tracker
//            TotalAbundance += ecosystemModelGrid.StateVariableGridTotal("Abundance", "NA", CohortFunctionalGroupIndices, cellIndices, "cohort", initialisation);
//
//            // Get total organic pool biomass
//            OrganicPoolOut = ecosystemModelGrid.GetEnviroGridTotal("Organic Pool", 0, cellIndices);
//
//            // Get total respiratory pool biomass
//            RespiratoryPoolOut = ecosystemModelGrid.GetEnviroGridTotal("Respiratory CO2 Pool", 0, cellIndices);
//
//            // Get total of all biomass
//            TotalBiomass = TotalLivingBiomass + RespiratoryPoolOut + OrganicPoolOut;
//
//            // Get number of cohorts and stocks
//            TotalNumberOfCohorts = globalDiagnosticVariables["NumberOfCohortsInModel"];
//            TotalNumberOfStocks = globalDiagnosticVariables["NumberOfStocksInModel"];
//
//            // Get numbers of cohort extinctions and productions
//            NumberOfCohortsExtinct = globalDiagnosticVariables["NumberOfCohortsExtinct"];
//            NumberOfCohortsProduced = globalDiagnosticVariables["NumberOfCohortsProduced"];
//            NumberOfCohortsCombined = globalDiagnosticVariables["NumberOfCohortsCombined"];
//        }
//
/** \brief
Generates the initial model outputs before the first time step
*/
@param ecosystemModelGrid The model grid 
@param cohortFunctionalGroupDefinitions The functional group definitions of cohorts in the model 
@param stockFunctionalGroupDefinitions The functional group definitions of stocks in the model 
@param cellIndices The list of indices of active cells in the model grid 
@param globalDiagnosticVariables A list of global diagnostic variables 
//         void InitialOutputs(ModelGrid ecosystemModelGrid, FunctionalGroupDefinitions cohortFunctionalGroupDefinitions, FunctionalGroupDefinitions
//            stockFunctionalGroupDefinitions, List<unsigned[]> cellIndices, map<string,double> globalDiagnosticVariables, MadingleyModelInitialisation initialisation)
//        {
//            // Calculate the output variables
//            CalculateOutputs(cohortFunctionalGroupDefinitions, stockFunctionalGroupDefinitions, ecosystemModelGrid, cellIndices, 
//                globalDiagnosticVariables, initialisation);
//
//            // Generate the initial console outputs
//            InitialConsoleOutputs();
//
//            // Generate the initial file outputs
//            InitialFileOutptus(ecosystemModelGrid);
//        }
//
/** \brief
Generates the initial console outputs
*/
//        private void InitialConsoleOutputs()
//        {
//            // Console outputs for all levels of detail
//            Console.ForegroundColor = ConsoleColor.Green;
//            Console.WriteLine("Initial living biomass = " + String.Format("{0:E}", TotalLivingBiomass / 1000) + " kg");
//            Console.WriteLine("Initial total abundance  = " + String.Format("{0:E}", TotalAbundance) + " inds / km^2");
//            // Console outputs for medium or high detail levels
//            if ((ModelOutputDetail == OutputDetailLevel.Medium) || (ModelOutputDetail == OutputDetailLevel.High))
//            {
//                Console.WriteLine("Initial total number of cohorts = " + TotalNumberOfCohorts);
//                Console.WriteLine("Initial number of stocks = " + TotalNumberOfStocks);
//
//                // Console outputs for high detail level
//                if (ModelOutputDetail == OutputDetailLevel.High)
//                {
//                    Console.WriteLine("Initial total biomass (all) = " + String.Format("{0:E}", TotalBiomass / 1000) + " kg", "high", ConsoleColor.White);
//                    Console.WriteLine("Initial respiratory pool biomass = " + String.Format("{0:E}", RespiratoryPoolOut / 1000) + " kg", "high", ConsoleColor.White);
//                    Console.WriteLine("Initial organic pool biomass = " + String.Format("{0:E}", OrganicPoolOut / 1000) + " kg", "high", ConsoleColor.White);
//                }
//            }
//            // Blank line at end of console initial outputs
//            Console.WriteLine(" ");
//            Console.ForegroundColor = ConsoleColor.White;
//        }
//
/** \brief
Generates the initial file outputs
*/
@param ecosystemModelGrid The model grid 
//        private void InitialFileOutptus(ModelGrid ecosystemModelGrid)
//        {
//            Console.WriteLine("Writing initial global outputs to file...");
//
//            // Write out total biomass in the organic and respiratory pools to the relevant one-dimensional output variables
//            DataConverter.ValueToSDS1D(OrganicPoolOut, "Organic matter pool", "Time step", ecosystemModelGrid.GlobalMissingValue,
//                BasicOutput, 0);
//            DataConverter.ValueToSDS1D(RespiratoryPoolOut, "Respiratory CO2 pool", "Time step", ecosystemModelGrid.GlobalMissingValue,
//                BasicOutput, 0);
//
//            // Write out the total number of cohorts and stocks to the relevant one-dimensional output variables
//            DataConverter.ValueToSDS1D(TotalNumberOfCohorts, "Number of cohorts in model", "Time step", ecosystemModelGrid.
//                GlobalMissingValue, BasicOutput, 0);
//            DataConverter.ValueToSDS1D(TotalNumberOfStocks, "Number of stocks in model", "Time step", ecosystemModelGrid.
//                GlobalMissingValue, BasicOutput, 0);
//
//            //Write out the number of cohorts produced, extinct and combined to the relevant 1D output variables
//            DataConverter.ValueToSDS1D(NumberOfCohortsProduced, "Number of cohorts produced", "Time step", ecosystemModelGrid.
//                GlobalMissingValue, BasicOutput, 0);
//            DataConverter.ValueToSDS1D(NumberOfCohortsExtinct, "Number of cohorts extinct", "Time step", ecosystemModelGrid.
//                GlobalMissingValue, BasicOutput, 0);
//            DataConverter.ValueToSDS1D(NumberOfCohortsCombined, "Number of cohorts combined", "Time step", ecosystemModelGrid.
//                GlobalMissingValue, BasicOutput, 0);
//
//            // Write out the total living biomass
//            DataConverter.ValueToSDS1D(TotalLivingBiomass, "Total living biomass", "Time step", ecosystemModelGrid.GlobalMissingValue,
//                BasicOutput, 0);
//        }
//
/** \brief
Generates the outputs for the current time step
*/
@param ecosystemModelGrid The model grid 
@param currentTimeStep The current model time step 
@param timeStepTimer The timer for the current time step 
@param cohortFunctionalGroupDefinitions The functional group definitions of cohorts in the model 
@param stockFunctionalGroupDefinitions The functional group definitions of stocks in the model 
@param cellIndices List of indices of active cells in the model grid 
@param globalDiagnosticVariables The global diagnostic variables for the model run 
//         void TimeStepOutputs(ModelGrid ecosystemModelGrid, unsigned currentTimeStep, unsigned currentMonth, StopWatch timeStepTimer, FunctionalGroupDefinitions
//            cohortFunctionalGroupDefinitions, FunctionalGroupDefinitions stockFunctionalGroupDefinitions, List<unsigned[]> cellIndices,
//            map<string,double> globalDiagnosticVariables, MadingleyModelInitialisation initialisation)
//        {
//            // Calculate the output variables for this time step
//            CalculateOutputs(cohortFunctionalGroupDefinitions, stockFunctionalGroupDefinitions, ecosystemModelGrid, cellIndices,
//                globalDiagnosticVariables, initialisation);
//
//            // Generate the time step console outputs
//            TimeStepConsoleOutputs(currentTimeStep, currentMonth, timeStepTimer);
//
//            // Generate the time step file outputs
//            TimeStepFileOutputs(ecosystemModelGrid, currentTimeStep);
//        }
//
/** \brief
Generates the console outputs for the current time step
*/
@param currentTimeStep The current model time step 
@param timeStepTimer The timer for the current time step 
//        private void TimeStepConsoleOutputs(unsigned currentTimeStep, unsigned currentMonth, StopWatch timeStepTimer)
//        {
//            // Console outputs for all levels of detail
//            Console.ForegroundColor = ConsoleColor.Red;
//            Console.WriteLine("Completed time step " + (currentTimeStep + 1) + "(Month: " + (currentMonth + 1) + ")\n");
//            Console.ForegroundColor = ConsoleColor.White;
//            Console.WriteLine("Elapsed time in seconds this time step: " + timeStepTimer.GetElapsedTimeSecs() + "\n");
//            Console.ForegroundColor = ConsoleColor.Green;
//            Console.WriteLine("Living biomass = " + String.Format("{0:E}", TotalLivingBiomass / 1000) + " kg");
//            Console.WriteLine("Total abundance  = " + String.Format("{0:E}", TotalAbundance) + " inds / km^2");
//            // Console outputs for medium or high detail levels
//            if ((ModelOutputDetail == OutputDetailLevel.Medium) || (ModelOutputDetail == OutputDetailLevel.High))
//            {
//                Console.WriteLine("Total number of cohorts = " + TotalNumberOfCohorts);
//                Console.WriteLine("Total number of stocks = " + TotalNumberOfStocks);
//                Console.WriteLine("Number of cohorts extinct = " + NumberOfCohortsExtinct);
//                Console.WriteLine("Number of cohorts produced = " + NumberOfCohortsProduced);
//                Console.WriteLine("Number of cohorts combined = " + NumberOfCohortsCombined);
//
//                // Console ouputs for high detail level
//                if (ModelOutputDetail == OutputDetailLevel.High)
//                {
//                    Console.WriteLine("Total biomass (all) = " + String.Format("{0:E}", TotalBiomass / 1000) + " kg");
//                    Console.WriteLine("Respiratory pool biomass = " + String.Format("{0:E}", RespiratoryPoolOut / 1000) + " kg");
//                    Console.WriteLine("Organic pool biomass = " + String.Format("{0:E}", OrganicPoolOut / 1000) + " kg");
//
//                }
//            }
//            // Blank line at end of console time step outputs
//            Console.WriteLine(" ");
//            Console.ForegroundColor = ConsoleColor.White;
//        }
//
/** \brief
Generates the file outputs for the current time step
*/
@param ecosystemModelGrid The model grid 
@param currentTimeStep The current model time step 
//        private void TimeStepFileOutputs(ModelGrid ecosystemModelGrid, unsigned currentTimeStep)
//        {
//            Console.WriteLine("Writing global ouputs to file...\n");
//
//            // Write out the basic diagnostic variables
//            DataConverter.ValueToSDS1D(OrganicPoolOut, "Organic matter pool", "Time step",
//                ecosystemModelGrid.GlobalMissingValue, BasicOutput, (int)currentTimeStep + 1);
//            DataConverter.ValueToSDS1D(RespiratoryPoolOut, "Respiratory CO2 pool", "Time step",
//                ecosystemModelGrid.GlobalMissingValue, BasicOutput, (int)currentTimeStep + 1);
//            DataConverter.ValueToSDS1D(TotalNumberOfCohorts, "Number of cohorts in model", "Time step",
//                ecosystemModelGrid.GlobalMissingValue, BasicOutput, (int)currentTimeStep + 1);
//            DataConverter.ValueToSDS1D(TotalNumberOfStocks, "Number of stocks in model", "Time step",
//                ecosystemModelGrid.GlobalMissingValue, BasicOutput, (int)currentTimeStep + 1);
//            DataConverter.ValueToSDS1D(NumberOfCohortsExtinct, "Number of cohorts extinct", "Time step",
//                ecosystemModelGrid.GlobalMissingValue, BasicOutput, (int)currentTimeStep + 1);
//            DataConverter.ValueToSDS1D(NumberOfCohortsProduced, "Number of cohorts produced", "Time step",
//                ecosystemModelGrid.GlobalMissingValue, BasicOutput, (int)currentTimeStep + 1);
//            DataConverter.ValueToSDS1D(NumberOfCohortsCombined, "Number of cohorts combined", "Time step",
//                ecosystemModelGrid.GlobalMissingValue, BasicOutput, (int)currentTimeStep + 1);
//            DataConverter.ValueToSDS1D(TotalLivingBiomass, "Total living biomass", "Time step",
//                ecosystemModelGrid.GlobalMissingValue, BasicOutput, (int)currentTimeStep + 1);
//        }
//
/** \brief
 Write to the output file values of the output variables at the end of the model run
*/
//         void FinalOutputs()
//        {
//            // Write out final total biomass to the console
//            Console.WriteLine("Final total biomass = {0}", TotalBiomass);
//
//        }
//
//    }
//}
#endif
