#ifndef MADINGLEYMODELINITIALISATION_H
#define MADINGLEYMODELINITIALISATION_H

#include <iostream>
#include<map>
#include <string>
#include <Properties.h>
#include <UtilityFunctions.h>
#include <FunctionalGroupDefinitions.h>
#include <EnviroData.h>
#include <MassBinsHandler.h>
/**
\file MadingleyModelInitialisation.h
\brief The MadingleyModelInitialisation header file
*/
//Changes:
//MB 7/8/2015
//replaced SortedList<> with map<>

//
//namespace Madingley
//{
/**
\brief
Initialization information for Madingley model simulations
*/
class MadingleyModelInitialisation
    {
    public:
/** String identifying time step units to be used by the simulations */

        StringProperty GlobalModelTimeStepUnit;
/** The number of time steps to be run in the simulations */
        UnsignedProperty NumTimeSteps;
/** The size of cells to be used in the model grid */
        FloatProperty CellSize;
/** The lowest extent of the model grid in degrees */
        FloatProperty BottomLatitude;
/** The uppermost extent of the model grid in degrees */
        FloatProperty TopLatitude;
/** The leftmost extent of the model grid in degrees */
        FloatProperty LeftmostLongitude;
/** The rightmost extent of the model grid in degrees */
        FloatProperty RightmostLongitude;
/** The rarefaction of active cells in the model grid */
        IntProperty CellRarefaction;
/** Whether to run the model for different grid cells in parallel (??) */
         BoolProperty RunInParallel;
/** Whether to run the model for different grid cells in parallel */
         BoolProperty RunCellsInParallel;
/** Whether to run the model for different simulations in parallel */
         BoolProperty RunSimulationsInParallel;
/** What's this? I dunno MB 6/8/2015 */
         StringProperty RunRealm;
/** Whether to draw cohort properties randomly when seeding them and whether cohorts will undergo ecological processes in a random order
\remark Value should be set in initialization file, but default value is true */
         BoolProperty DrawRandomly;
         //MB 6/8/2015 Now set in the constructor
/** The threshold abundance below which cohorts will be made extinct */
         DoubleProperty ExtinctionThreshold;
/** The threshold difference between cohorts, within which they will be merged */
         DoubleProperty MergeDifference;
/** The maximum number of cohorts to be in the model, per grid cell, when it is running */
         IntProperty MaxNumberOfCohorts;
/**Whether to run only dispersal (i.e. turn all other ecological processes off, and set dispersal probability to one temporarily) */
         BoolProperty DispersalOnly;
         //MB 6/8/2015 Default now set in the constructor
/** The weight threshold (grams) below which marine organisms that are not obligate zooplankton will be dispersed planktonically */
         DoubleProperty PlanktonDispersalThreshold;
/** The full path for the output files for a set of simulations */
           StringProperty OutputPath;
/** Whether to output detailed diagnostics for the ecological processes */
           BoolProperty TrackProcesses;
           //MB 6/8/2015 default set in constructor
/** Whether to output detailed diagnostics for the ecological processes*/
           BoolProperty TrackGlobalProcesses;
           //MB 6/8/2015 default set in constructor
/** Whether to display live outputs using Dataset Viewer during the model runs */
           BoolProperty LiveOutputs;
/** Whether or not to track trophic level biomass and flow information specific to the marine realm */
           BoolProperty TrackMarineSpecifics;

/** \brief Information from the initialization file
*/
        map<string, string> InitialisationFileStrings;
//        {
//            get { return _InitialisationFileStrings; }
//            set { _InitialisationFileStrings = value; }
//        }
/** \brief The functional group definitions of cohorts in the model */
          FunctionalGroupDefinitions _CohortFunctionalGroupDefinitions;

//        {
//            get { return _CohortFunctionalGroupDefinitions; }
//            set { _CohortFunctionalGroupDefinitions = value; }
//        }
/** \brief The functional group definitions of stocks in the model */
          FunctionalGroupDefinitions StockFunctionalGroupDefinitions;

//        {
//            get { return _StockFunctionalGroupDefinitions; }
//            set { _StockFunctionalGroupDefinitions = value; }
//        }
/** \brief The environmental layers for use in the model */
        map<string, EnviroData> EnviroStack; 
//        {
//            get { return _EnviroStack; }
//            set { _EnviroStack = value; }
//        }
/** \brief The paths and filenames for the diagnostics for the ecological processes */

          map<string, string> ProcessTrackingOutputs;
//        {
//            get { return _ProcessTrackingOutputs; }
//            set { _ProcessTrackingOutputs = value; }
//        }
/** \brief The string values for the units of each environmental data layer */
          map<string,string> Units;
//        {
//            get { return _Units; }
//            set { _Units = value; }
//        }
/** \brief An instance of the mass bin handler for the current model run */
          MassBinsHandler ModelMassBins;
//        { get { return _ModelMassBins; } }
//
/** Instance of Utilities for timestep conversions */
          UtilityFunctions Utilities;
//
/**
\brief Reads the initalization file to get information for the set of simulations to be run
@param initialisationFile The name of the initialization file with information on the simulations to be run
@param outputPath The path to folder in which outputs will be stored */
        MadingleyModelInitialisation(string initialisationFile, string outputPath)
        {
//            // Write to console
            cout<<"Initializing model...\n"<<endl;
            //set default value of cohort random draw properties
            DrawRandomly=true;
            //set default dispersal
            DispersalOnly=false;
            //defaults for process tracking
            TrackProcesses=false;
            TrackGlobalProcesses=false;
//Class variable doesn't need initialising
//            // Initialize the mass bins to be used during the model run
//            _ModelMassBins = new MassBinsHandler();
//
            // Read the intialisation files and copy them to the output directory
            ReadAndCopyInitialisationFiles(initialisationFile, outputPath);
//
//            // Copy parameter values to an output file
            CopyParameterValues(outputPath);
//
//           
//
        }
//
/** \brief Reads in all initialisation files and copies them to the output directory for future reference 
@param initialisationFile The name of the initialization file with information on the simulations to be run
@param outputPath The path to folder in which outputs will be stored
//        /// <todo>Need to adjust this file to deal with incorrect inputs, extra columns etc by throwing an error</todo>
//        /// <todo>Also need to strip leading spaces</todo>*/
          void ReadAndCopyInitialisationFiles(string initialisationFile, string outputPath)
        {
            NumTimeSteps=5;
//            // Construct file name
//            string FileString = "msds:csv?file=input/Model setup/" + initialisationFile + "&openMode=readOnly";
//
//            // Copy the initialisation file to the output directory
//            System.IO.File.Copy("input/Model setup/" + initialisationFile, outputPath + initialisationFile, true);
//
//            // Read in the data
//            DataSet InternalData = DataSet.Open(FileString);
//
//            // Get the names of parameters in the initialization file
//            var VarParameters = InternalData.Variables[1].GetData();
//
//            // Get the values for the parameters
//            var VarValues = InternalData.Variables[0].GetData();
//
//            // Loop over the parameters
//            for (int row = 0; row < VarParameters.Length; row++)
//            {
//                // Switch based on the name of the parameter, and write the value to the appropriate field
//                switch (VarParameters.GetValue(row).ToString().ToLower())
//                {
//                    case "timestep units":
//                        _GlobalModelTimeStepUnit = VarValues.GetValue(row).ToString();
//                        break;
//                    case "length of simulation (years)":
//                        _NumTimeSteps = (uint)Utilities.ConvertTimeUnits("year",  _GlobalModelTimeStepUnit)*Convert.ToUInt32(VarValues.GetValue(row));
//                        break;
//                    case "number timesteps":
//                        _NumTimeSteps = Convert.ToUInt32(VarValues.GetValue(row));
//                        break;
//                    case "grid cell size":
//                        _CellSize = Convert.ToSingle(VarValues.GetValue(row));
//                        break;
//                    case "bottom latitude":
//                        _BottomLatitude = Convert.ToSingle(VarValues.GetValue(row));
//                        break;
//                    case "top latitude":
//                        _TopLatitude = Convert.ToSingle(VarValues.GetValue(row));
//                        break;
//                    case "leftmost latitude":
//                        _LeftmostLongitude = Convert.ToSingle(VarValues.GetValue(row));
//                        break;
//                    case "rightmost latitude":
//                        _RightmostLongitude = Convert.ToSingle(VarValues.GetValue(row));
//                        break;
//                    case "grid cell rarefaction":
//                        _CellRarefaction = Convert.ToInt32(VarValues.GetValue(row));
//                        Debug.Assert(_CellRarefaction >= 1, "Cell rarefaction cannot be less than 1");
//                        break;
//                    case "run cells in parallel":
//                        switch (VarValues.GetValue(row).ToString().ToLower())
//                        {
//                            case "yes":
//                                _RunCellsInParallel = true;
//                                break;
//                            case "no":
//                                _RunCellsInParallel = false;
//                                break;
//                        }
//                        break;
//                    case "run simulations in parallel":
//                        switch (VarValues.GetValue(row).ToString().ToLower())
//                        {
//                            case "yes":
//                                _RunSimulationsInParallel = true;
//                                break;
//                            case "no":
//                                _RunSimulationsInParallel = false;
//                                break;
//                        }
//                        break;
//                    case "run single realm":
//                        _RunRealm = VarValues.GetValue(row).ToString().ToLower();
//                        break;
//                    case "draw randomly":
//
//                        switch (VarValues.GetValue(row).ToString().ToLower())
//                        {
//                            case "yes":
//                                _DrawRandomly = true;
//                                break;
//                            case "no":
//                                _DrawRandomly = false;
//                                break;
//                        }
//                        break;
//                    case "extinction threshold":
//                        _ExtinctionThreshold = Convert.ToDouble(VarValues.GetValue(row));
//                        break;
//                    case "merge difference":
//                        _MergeDifference = Convert.ToDouble(VarValues.GetValue(row));
//                        break;
//                    case "maximum number of cohorts":
//                        _MaxNumberOfCohorts = Convert.ToInt32(VarValues.GetValue(row));
//                        break;
//                    case "track processes":
//                        switch (VarValues.GetValue(row).ToString().ToLower())
//                        {
//                            case "yes":
//                                _TrackProcesses = true;
//                                break;
//                            case "no":
//                                _TrackProcesses = false;
//                                break;
//                        }
//                        break;
//                    case "track global processes":
//                        switch (VarValues.GetValue(row).ToString().ToLower())
//                        {
//                            case "yes":
//                                _TrackGlobalProcesses = true;
//                                break;
//                            case "no":
//                                _TrackGlobalProcesses = false;
//                                break;
//                        }
//                        break;
//                    case "mass bin filename":
//                        // Set up the mass bins as specified in the initialization file
//                        _ModelMassBins.SetUpMassBins(VarValues.GetValue(row).ToString(), outputPath);
//                        break;
//                    case "new cohorts filename":
//                        _ProcessTrackingOutputs.Add("NewCohortsOutput", VarValues.GetValue(row).ToString());
//                        break;
//                    case "maturity filename":
//                        _ProcessTrackingOutputs.Add("MaturityOutput", VarValues.GetValue(row).ToString());
//                        break;
//                    case "biomasses eaten filename":
//                        _ProcessTrackingOutputs.Add("BiomassesEatenOutput", VarValues.GetValue(row).ToString());
//                        break;
//                    case "trophic flows filename":
//                        _ProcessTrackingOutputs.Add("TrophicFlowsOutput", VarValues.GetValue(row).ToString());
//                        break;
//                    case "growth filename":
//                        _ProcessTrackingOutputs.Add("GrowthOutput", VarValues.GetValue(row).ToString());
//                        break;
//                    case "metabolism filename":
//                        _ProcessTrackingOutputs.Add("MetabolismOutput", VarValues.GetValue(row).ToString());
//                        break;
//                    case "npp filename":
//                        _ProcessTrackingOutputs.Add("NPPOutput", VarValues.GetValue(row).ToString());
//                        break;
//                    case "predation flows filename":
//                        _ProcessTrackingOutputs.Add("PredationFlowsOutput", VarValues.GetValue(row).ToString());
//                        break;
//                    case "herbivory flows filename":
//                        _ProcessTrackingOutputs.Add("HerbivoryFlowsOutput",VarValues.GetValue(row).ToString());
//                        break;
//                    case "mortality filename":
//                        _ProcessTrackingOutputs.Add("MortalityOutput", VarValues.GetValue(row).ToString());
//                        break;
//                    case "extinction filename":
//                        _ProcessTrackingOutputs.Add("ExtinctionOutput", VarValues.GetValue(row).ToString());
//                        break;
//                    case "environmental data file":
//                        _InitialisationFileStrings.Add("Environmental", VarValues.GetValue(row).ToString());
//                        // Read environmental data layers
//                        this.ReadEnvironmentalLayers(VarValues.GetValue(row).ToString(), outputPath);
//                        break;
//                    case "cohort functional group definitions file":
//                        Console.WriteLine("Reading functional group definitions...\n");
//                        _InitialisationFileStrings.Add("CohortFunctional", VarValues.GetValue(row).ToString());
//                        // Open a the specified csv file and set up the cohort functional group definitions
//                        _CohortFunctionalGroupDefinitions = new FunctionalGroupDefinitions(VarValues.GetValue(row).ToString(), outputPath);
//                        break;
//                    case "stock functional group definitions file":
//                        _InitialisationFileStrings.Add("StockFunctional", VarValues.GetValue(row).ToString());
//                        // Open a the specified csv file and set up the stock functional group definitions
//                        _StockFunctionalGroupDefinitions = new FunctionalGroupDefinitions(VarValues.GetValue(row).ToString(), outputPath);
//                        break;
//                    case "specific location file":
//                        if (VarValues.GetValue(row).ToString() != "")
//                        {
//                            _InitialisationFileStrings.Add("Locations", VarValues.GetValue(row).ToString());
//                            // Copy the initialisation file to the output directory
//                            System.IO.File.Copy("input/Model setup/" + _InitialisationFileStrings["Locations"], outputPath + _InitialisationFileStrings["Locations"], true);
//                        }
//                        break;
//                    case "output detail":
//                        _InitialisationFileStrings.Add("OutputDetail", VarValues.GetValue(row).ToString());
//                        break;
//                    case "human npp extraction":
//                        _InitialisationFileStrings.Add("HumanNPPExtraction", VarValues.GetValue(row).ToString());
//                        break;
//                    case "dispersal only":
//                        if (VarValues.GetValue(row).ToString() == "yes")
//                            _DispersalOnly = true;
//                        else _DispersalOnly = false;
//                        break;
//                    case "dispersal only type":
//                        _InitialisationFileStrings.Add("DispersalOnlyType", VarValues.GetValue(row).ToString());
//                        break;
//                    case "plankton size threshold":
//                        _PlanktonDispersalThreshold =  Convert.ToDouble(VarValues.GetValue(row));
//                        break;
//                    case "live outputs":
//                        if (VarValues.GetValue(row).ToString() == "yes")
//                            _LiveOutputs = true;
//                        else _LiveOutputs = false;
//                        break;
//                    case "track marine specifics":
//                        if (VarValues.GetValue(row).ToString() == "yes")
//                            _TrackMarineSpecifics = true;
//                        else _TrackMarineSpecifics = false;
//                        break;
//                }
//
//            }
//
//            InternalData.Dispose();
        }
//
/** \brief Copy parameter values to a text file in the specified output directory
@param outputDirectory The directory for outputs */
void CopyParameterValues(string outputDirectory)
        {
//            // Create a stream write object to write the parameter values to
//            StreamWriter sw = new StreamWriter(outputDirectory + "Parameters.txt");
//   
//            // Write out the column headings
//            sw.WriteLine("Ecological process\tParameter name\tParameter value");
//
//            // Create dummy instances of the ecological processes
//            RevisedHerbivory DummyHerbivory = new RevisedHerbivory(0.0, _GlobalModelTimeStepUnit);
//            RevisedPredation DummyPredation = new RevisedPredation(0.0, _GlobalModelTimeStepUnit);
//            MetabolismEndotherm DummyEndoMetabolism = new MetabolismEndotherm(_GlobalModelTimeStepUnit);
//            MetabolismEctotherm DummyEctoMetabolism = new MetabolismEctotherm(_GlobalModelTimeStepUnit);
//            BackgroundMortality DummyBackgroundMortality = new BackgroundMortality(_GlobalModelTimeStepUnit);
//            SenescenceMortality DummySenescenceMortality = new SenescenceMortality(_GlobalModelTimeStepUnit);
//            StarvationMortality DummyStarvationMortality = new StarvationMortality(_GlobalModelTimeStepUnit);
//            ReproductionBasic DummyReproduction = new ReproductionBasic(_GlobalModelTimeStepUnit, _DrawRandomly);
//            DiffusiveDispersal DummyDiffusiveDispersal = new DiffusiveDispersal(_GlobalModelTimeStepUnit, _DrawRandomly);
//            RevisedTerrestrialPlantModel DummyPlantModel = new RevisedTerrestrialPlantModel();
//            Activity DummyActivityModel = new Activity();
//
//            
//            // Call the methods in these processes that write the parameter values out
//            DummyHerbivory.WriteOutParameterValues(sw);
//            DummyPredation.WriteOutParameterValues(sw);
//            DummyEndoMetabolism.WriteOutParameterValues(sw);
//            DummyEctoMetabolism.WriteOutParameterValues(sw);
//            DummyBackgroundMortality.WriteOutParameterValues(sw);
//            DummySenescenceMortality.WriteOutParameterValues(sw);
//            DummyStarvationMortality.WriteOutParameterValues(sw);
//            DummyReproduction.WriteOutParameterValues(sw);
//            DummyDiffusiveDispersal.WriteOutParameterValues(sw);
//            DummyPlantModel.WriteOutParameterValues(sw);
//            DummyActivityModel.WriteOutParameterValues(sw);
//
//
//            sw.Dispose();
//
        }
//

/** \brief Reads the environmental layers listed in the specified file containing a list of environmental layers
//        /// </summary>
@param environmentalLayerFile The name of the file containing the list of environmental layers
@param outputPath The path to folder in which outputs will be stored */
void ReadEnvironmentalLayers(string environmentalLayerFile, string outputPath){

            cout<<"Reading in environmental data:"<<endl;
//
//            // Declare lists to hold the information required to read the environmental layers
//            List<string> Folders = new List<string>();
//            List<string> Filenames = new List<string>();
//            List<string> DatasetNames = new List<string>();
//            List<string> FileTypes = new List<string>();
//            List<string> LayerName = new List<string>();
//            List<string> StaticLayer = new List<string>();
//            List<string> Extensions = new List<string>();
//            List<string> Resolutions = new List<string>();
//            List<string> MethodUnits = new List<string>();
//
//            // Variable to store the file name of the environmental data files
//            string TempFilename;
//
//            // Construct the full URI for the file  containing the list of environmental layers
//            string FileString = "msds:csv?file=input/Model setup/" + environmentalLayerFile + "&openMode=readOnly";
//
//            //Copy the file containing the list of environmental layers to the output directory
//            System.IO.File.Copy("input/Model setup/" + environmentalLayerFile, outputPath + environmentalLayerFile, true);
//
//            // Read in the data
//            DataSet InternalData = DataSet.Open(FileString);
//
//            // Loop over the parameters associated with the list of environmental layers
//           
//            foreach (Variable v in InternalData.Variables)
//            {
//                // Get the name of the parameter
//                string HeaderName = v.Name;
//
//                // Create a local copy of all of the values associated with this parameter
//                var TempValues = v.GetData();
//
//                // Switch based on the name of the parameter, and store the parameter values in the appropriate list
//                switch (HeaderName.ToLower())
//                {
//                    case "folder":
//                        for (int ii = 0; ii < TempValues.Length; ii++) Folders.Add(TempValues.GetValue(ii).ToString());
//                        break;
//                    case "filename":
//                        for (int ii = 0; ii < TempValues.Length; ii++) Filenames.Add(TempValues.GetValue(ii).ToString());
//                        break;
//                    case "extension":
//                        for (int ii = 0; ii < TempValues.Length; ii++) Extensions.Add(TempValues.GetValue(ii).ToString());
//                        break;
//                    case "dataset name":
//                        for (int ii = 0; ii < TempValues.Length; ii++) DatasetNames.Add(TempValues.GetValue(ii).ToString());
//                        break;
//                    case "filetype":
//                        for (int ii = 0; ii < TempValues.Length; ii++) FileTypes.Add(TempValues.GetValue(ii).ToString());
//                        break;
//                    case "internal layer name":
//                        for (int ii = 0; ii < TempValues.Length; ii++) LayerName.Add(TempValues.GetValue(ii).ToString());
//                        break;
//                    case "static":
//                        for (int ii = 0; ii < TempValues.Length; ii++) StaticLayer.Add(TempValues.GetValue(ii).ToString());
//                        break;
//                    case "resolution":
//                        for (int ii = 0; ii < TempValues.Length; ii++) Resolutions.Add(TempValues.GetValue(ii).ToString());
//                        break;
//                    case "units":
//                        for (int ii = 0; ii < TempValues.Length; ii++)
//                        {
//                            MethodUnits.Add(TempValues.GetValue(ii).ToString());
//                        }
//                        break;
//                }
//            }
//
//            for (int ii = 0; ii < MethodUnits.Count; ii++)
//            {
//                Units.Add(LayerName[ii], MethodUnits[ii]);
//            }
//
//            // Check that there are the same number of values for all parameters
//            Debug.Assert(Folders.Count() == Filenames.Count() && Filenames.Count() == DatasetNames.Count() && DatasetNames.Count() == FileTypes.Count() && FileTypes.Count() == LayerName.Count(),
//                "Error in Environmental Data Layer import lists - unequal number of filenames, dataset names, filetypes and datalayer names");
//
//            // Loop over parameter values
//            for (int ii = 0; ii < Filenames.Count(); ii++)
//            {
//                Console.Write("\rVariable {0} of {1} ", ii+1,Filenames.Count);
//                // If the layers are not static, then suffix the file name with '1' - not currently implemented
//                if (StaticLayer[ii].ToLower().Equals("n"))
//                {
//                    Debug.Fail("This option is currently not supported");
//                    Filenames[ii] = Filenames[ii] + "1";
//                }
//                // For layers where the file format is ESRI ASCII grid, the dataset name is the same as the file name
//                if (FileTypes[ii].ToLower().Equals("esriasciigrid"))
//                {
//                    DatasetNames[ii] = Filenames[ii];
//                }
//                // Generate the appropriate file name for the environmental data layer
//                if (Folders[ii].ToLower().Equals("input"))
//                {
//                    TempFilename = "input/Data/" + Filenames[ii];
//                }
//                else
//                {
//                    TempFilename = Folders[ii] + "/" + Filenames[ii];
//                }
//                Filenames[ii] = TempFilename + Extensions[ii];
//                // Read in and store the environmental data
//                EnviroStack.Add(LayerName[ii], new EnviroData(Filenames[ii], DatasetNames[ii], FileTypes[ii], Resolutions[ii], MethodUnits[ii]));
//            }
//            Console.WriteLine("\n\n");
//        }
//
//
    }
};
#endif