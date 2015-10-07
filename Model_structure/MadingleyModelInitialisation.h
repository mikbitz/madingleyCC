#ifndef MADINGLEYMODELINITIALISATION_H
#define MADINGLEYMODELINITIALISATION_H

#include <iostream>
#include <map>
#include <string>
#include <UtilityFunctions.h>
#include <FunctionalGroupDefinitions.h>
#include <EnviroData.h>
#include <MassBinsHandler.h>
#include <sstream>
/**
 \ file *MadingleyModelInitialisation.h
 \brief The MadingleyModelInitialisation header file
 */

/**
 \ brief*
 Initialization information for Madingley model simulations
 */
class MadingleyModelInitialisation {
public:
    //----------------------------------------------------------------------------------------------
    //Variables
    //----------------------------------------------------------------------------------------------
     /** String identifying time step units to be used by the simulations */
    string GlobalModelTimeStepUnit;
    /** The number of time steps to be run in the simulations */
    unsigned NumTimeSteps;
    /** The size of cells to be used in the model grid */
    float CellSize;
    /** The lowest extent of the model grid in degrees */
    float BottomLatitude;
    /** The uppermost extent of the model grid in degrees */
    float TopLatitude;
    /** The leftmost extent of the model grid in degrees */
    float LeftmostLongitude;
    /** The rightmost extent of the model grid in degrees */
    float RightmostLongitude;
    /** The rarefaction of active cells in the model grid */
    int CellRarefaction;
    /** Whether to run the model for different grid cells in parallel (??) */
    bool RunInParallel;
    /** Whether to run the model for different grid cells in parallel */
    bool RunCellsInParallel;
    /** Whether to run the model for different simulations in parallel */
    bool RunSimulationsInParallel;
    /** Include marine realm ? */
    string RunRealm;
    /** Whether to draw cohort properties randomly when seeding them and whether cohorts will undergo ecological processes in a random order
     @remark Value should be set in initialization file, but default value is true */
    bool DrawRandomly;
    //MB 6/8/2015 Now set in the constructor
    /** The threshold abundance below which cohorts will be made extinct */
    double ExtinctionThreshold;
    /** The threshold difference between cohorts, within which they will be merged */
    double MergeDifference;
    /** The maximum number of cohorts to be in the model, per grid cell, when it is running */
    int MaxNumberOfCohorts;
    /**Whether to run only dispersal (i.e. turn all other ecological processes off, and set dispersal probability to one temporarily) */
    bool DispersalOnly;
    //MB 6/8/2015 Default now set in the constructor
    /** The weight threshold (grams) below which marine organisms that are not obligate zooplankton will be dispersed planktonically */
    double PlanktonDispersalThreshold;
    /** The full path for the output files for a set of simulations */
    string OutputPath;
    /** Whether to output detailed diagnostics for the ecological processes */
    bool TrackProcesses;
    //MB 6/8/2015 default set in constructor
    /** Whether to output detailed diagnostics for the ecological processes*/
    bool TrackGlobalProcesses;
    //MB 6/8/2015 default set in constructor
    /** Whether to display live outputs using Dataset Viewer during the model runs */
    bool LiveOutputs;
    /** Whether or not to track trophic level biomass and flow information specific to the marine realm */
    bool TrackMarineSpecifics;
    /** \brief Information from the initialization file  */
    map<string, string> InitialisationFileStrings;
    /** \brief The functional group definitions of cohorts in the model */
    FunctionalGroupDefinitions CohortFunctionalGroupDefinitions;
    /** \brief The functional group definitions of stocks in the model */
    FunctionalGroupDefinitions StockFunctionalGroupDefinitions;
    /** \brief The environmental layers for use in the model */
    map<string, EnviroData> EnviroStack;
    /** \brief The paths and filenames for the diagnostics for the ecological processes */
    map<string, string> ProcessTrackingOutputs;
    /** \brief The string values for the units of each environmental data layer */
    map<string, string> Units;
    /** \brief An instance of the mass bin handler for the current model run */
    MassBinsHandler ModelMassBins;
    /** Instance of Utilities for timestep conversions */
    UtilityFunctions Utilities;

    //----------------------------------------------------------------------------------------------
    //Methods
    //----------------------------------------------------------------------------------------------
 
    //----------------------------------------------------------------------------------------------
    /** \brief Reads the initialization file to get information for the set of simulations to be run
     @param initialisationFile The name of the initialization file with information on the simulations to be run
     @param outputPath The path to folder in which outputs will be stored */
    MadingleyModelInitialisation(string initialisationFile, string outputPath) {
        //            // Write to console
        cout << "Initializing model...\n" << endl;
        //set default value of cohort random draw properties
        DrawRandomly = true;
        //set default dispersal
        DispersalOnly = false;
        //defaults for process tracking
        TrackProcesses = false;
        TrackGlobalProcesses = false;

        // Read the intialisation files and copy them to the output directory
        ReadAndCopyInitialisationFiles(initialisationFile, outputPath);
    }
    //----------------------------------------------------------------------------------------------
    /** \brief Reads in all initialisation files and copies them to the output directory for future reference 
     @ param* initialisationFile The name of the initialization file with information on the simulations to be run
     @param outputPath The path to folder in which outputs will be stored
     //        /// <todo>Need to adjust this file to deal with incorrect inputs, extra columns etc by throwing an error</todo>
     //        /// <todo>Also need to strip leading spaces</todo>*/
    void ReadAndCopyInitialisationFiles(string initialisationFile, string outputPath) {
        NumTimeSteps = 5;
        cout << "Reading initialisation parameters file..." << endl;
        ifstream infile(initialisationFile.c_str());
        if (infile.is_open()) {
            string l, header[2];
            getline(infile, l);
            //trim off newline character
            l.pop_back();
            istringstream s(l);
            for (unsigned i = 0; i < 2; i++) {
                getline(s, header[i], ',');
                transform(header[i].begin(), header[i].end(), header[i].begin(), ::tolower);
            }
            while (infile.good()) {
                string l, data[2];
                getline(infile, l);
                //trim off newline if there is still more data
                if (infile.good())l.pop_back();
                if (l.length() > 1) {
                    istringstream s(l);
                    for (unsigned i = 0; i < 2; i++) {
                        getline(s, data[i], ',');
                    }
                }
                string param = data[0], val = data[1];
                transform(param.begin(), param.end(), param.begin(), ::tolower);
                transform(val.begin(), val.end(), val.begin(), ::tolower);

                if (param == "timestep units")GlobalModelTimeStepUnit = val;
                if (param == "length of simulation (years)")NumTimeSteps = (unsigned) Utilities.ConvertTimeUnits("year", GlobalModelTimeStepUnit) * atoi(val.c_str());
                if (param == "number timesteps") NumTimeSteps = atoi(val.c_str());
                if (param == "grid cell size") CellSize = atof(val.c_str());
                if (param == "bottom latitude") BottomLatitude = atof(val.c_str());
                if (param == "top latitude") TopLatitude = atof(val.c_str());
                if (param == "leftmost longitude") LeftmostLongitude = atof(val.c_str());
                if (param == "rightmost longitude") RightmostLongitude = atof(val.c_str());
                if (param == "grid cell rarefaction") CellRarefaction = atoi(val.c_str());
                if (param == "run cells in parallel") RunCellsInParallel = ((val == "yes") ? true : false);
                if (param == "run simulations in parallel") RunSimulationsInParallel = ((val == "yes") ? true : false);
                if (param == "run single realm") RunRealm = val;
                if (param == "draw randomly") DrawRandomly = ((val == "yes") ? true : false);
                if (param == "extinction threshold") ExtinctionThreshold = atof(val.c_str());
                if (param == "merge difference") MergeDifference = atof(val.c_str());
                if (param == "maximum number of cohorts") MaxNumberOfCohorts = atof(val.c_str());
                if (param == "dispersal only") DispersalOnly = ((val == "yes") ? true : false);
                if (param == "plankton size threshold") PlanktonDispersalThreshold = atof(val.c_str());
                if (param == "live outputs") LiveOutputs = ((val == "yes") ? true : false);
                if (param == "track marine specifics") TrackMarineSpecifics = ((val == "yes") ? true : false);
                if (param == "track processes") TrackProcesses = ((val == "yes") ? true : false);
                if (param == "track global processes") TrackGlobalProcesses = ((val == "yes") ? true : false);

                if (param == "new cohorts filename") ProcessTrackingOutputs["NewCohortsOutput"] = val;
                if (param == "maturity filename") ProcessTrackingOutputs["MaturityOutput"] = val;
                if (param == "biomasses eaten filename") ProcessTrackingOutputs["BiomassesEatenOutput"] = val;
                if (param == "trophic flows filename") ProcessTrackingOutputs["TrophicFlowsOutput"] = val;
                if (param == "growth filename") ProcessTrackingOutputs["GrowthOutput"] = val;
                if (param == "metabolism filename") ProcessTrackingOutputs["MetabolismOutput"] = val;
                if (param == "npp filename") ProcessTrackingOutputs["NPPOutput"] = val;
                if (param == "predation flows filename") ProcessTrackingOutputs["PredationFlowsOutput"] = val;
                if (param == "herbivory flows filename") ProcessTrackingOutputs["HerbivoryFlowsOutput"] = val;
                if (param == "mortality filename") ProcessTrackingOutputs["MortalityOutput"] = val;
                if (param == "extinction filename") ProcessTrackingOutputs["ExtinctionOutput"] = val;

                if (param == "output detail") InitialisationFileStrings["OutputDetail"] = val;
                if (param == "human npp extraction") InitialisationFileStrings["HumanNPPExtraction"] = val;
                if (param == "dispersal only type") InitialisationFileStrings["DispersalOnlyType"] = val;
                if (param == "cohort functional group definitions file") InitialisationFileStrings["CohortFunctional"] = val;
                if (param == "stock functional group definitions file") InitialisationFileStrings["StockFunctional"] = val;
                if (param == "specific location file") InitialisationFileStrings["Locations"] = val;
                if (param == "environmental data file") InitialisationFileStrings["Environmental"] = val;

                //read in mass bins : use data[1] so that filename isn't lower-cased
                if (param == "mass bin filename") ModelMassBins.SetUpMassBins(data[1], outputPath);

                //Read environmental data layers
                if (param == "environmental data file") {
                    ReadEnvironmentalLayers(data[1], outputPath);
                }
                if (param == "cohort functional group definitions file") {
                    cout << "Reading functional group definitions..." << endl;
                    CohortFunctionalGroupDefinitions = FunctionalGroupDefinitions(data[1], outputPath);
                }
                if (param == "stock functional group definitions file") {
                    cout << "Reading stock group definitions..." << endl;

                    //Open a the specified csv file and set up the stock functional group definitions
                    StockFunctionalGroupDefinitions = FunctionalGroupDefinitions(data[1], outputPath);
                }



            }
        } else {
            cout << "Something wrong with initialisation parameter file " << initialisationFile << endl;
        }

        assert(CellRarefaction >= 1 && "Cell rarefaction cannot be less than 1");
    }
    //----------------------------------------------------------------------------------------------
    /** \brief Reads the environmental layers listed in the specified file containing a list of environmental layers
     @param environmentalLayerFile The name of the file containing the list of environmental layers
     @param outputPath The path to folder in which outputs will be stored */
    void ReadEnvironmentalLayers(string environmentalLayerFile, string outputPath) {

        cout << "Reading in environmental data:" << endl;
        ifstream infile(environmentalLayerFile.c_str());
        //            // Declare lists to hold the information required to read the environmental layers
        vector<string> Folders;
        vector<string> Filenames;
        vector<string> DatasetNames;
        vector<string> FileTypes;
        vector<string> LayerName;
        vector<string> StaticLayer;
        vector<string> Extensions;
        vector<string> Resolutions;
        vector<string> MethodUnits;
        if (infile.is_open()) {

            string l, header[9];
            getline(infile, l);
            //trim off newline character
            l.pop_back();
            istringstream s(l);
            for (unsigned i = 0; i < 9; i++) {
                getline(s, header[i], ',');
                transform(header[i].begin(), header[i].end(), header[i].begin(), ::tolower);
                if (header[i] != "folder" &&
                        header[i] != "filename" &&
                        header[i] != "extension" &&
                        header[i] != "dataset name" &&
                        header[i] != "filetype" &&
                        header[i] != "internal layer name" &&
                        header[i] != "static" &&
                        header[i] != "resolution" &&
                        header[i] != "units") {
                    cout << "Bad header in environmentalLayerFile file " << environmentalLayerFile << endl;
                    exit(1);
                }

            }


            while (infile.good()) {
                string l, data;
                getline(infile, l);
                if (infile.good())l.pop_back();
                if (l.length() > 1) {
                    istringstream s(l);
                    for (unsigned i = 0; i < 9; i++) {
                        getline(s, data, ',');

                        if (header[i] == "folder") Folders.push_back(data);
                        if (header[i] == "filename") Filenames.push_back(data);
                        if (header[i] == "extension") Extensions.push_back(data);
                        if (header[i] == "dataset name") DatasetNames.push_back(data);
                        if (header[i] == "filetype") FileTypes.push_back(data);
                        if (header[i] == "internal layer name") LayerName.push_back(data);
                        if (header[i] == "static") StaticLayer.push_back(data);
                        if (header[i] == "resolution") Resolutions.push_back(data);
                        if (header[i] == "units") MethodUnits.push_back(data);

                    }
                }
            }

        } else {
            cout << "Something wrong with environment parameter file " << environmentalLayerFile << endl;
            exit(1);
        }
        infile.close();


        for (int ii = 0; ii < MethodUnits.size(); ii++) {
            Units[LayerName[ii]] = MethodUnits[ii];
        }

        // Check that there are the same number of values for all parameters
        assert(Folders.size() == Filenames.size() && Filenames.size() == DatasetNames.size() && DatasetNames.size() == FileTypes.size() && FileTypes.size() == LayerName.size() &&
                "Error in Environmental Data Layer import lists - unequal number of filenames, dataset names, filetypes and datalayer names");

        // Loop over parameter values
        for (int ii = 0; ii < Filenames.size(); ii++) {
            cout << Filenames[ii] << ": Variable " << ii + 1 << "  of " << Filenames.size() << endl;
            // If the layers are not static, then suffix the file name with '1' - not currently implemented
            if (StaticLayer[ii] == "n") {
                assert("This option is currently not supported");
                Filenames[ii] = Filenames[ii] + "1";
            }
            // For layers where the file format is ESRI ASCII grid, the dataset name is the same as the file name
            if (FileTypes[ii] == "esriasciigrid") {
                DatasetNames[ii] = Filenames[ii];
            }
            // Generate the appropriate file name for the environmental data layer
            if (Folders[ii] == "input") {
                Filenames[ii] = "input/Data/" + Filenames[ii];
            } else {
                Filenames[ii] = Folders[ii] + "/" + Filenames[ii];
            }
            Filenames[ii] = Filenames[ii] + Extensions[ii];
            // Read in and store the environmental data
            EnviroStack[LayerName[ii]] = EnviroData(Filenames[ii], DatasetNames[ii], FileTypes[ii], Resolutions[ii], MethodUnits[ii]);
        }
        cout << endl;

    }
    //----------------------------------------------------------------------------------------------    
};
#endif