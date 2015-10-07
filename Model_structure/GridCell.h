#ifndef GRIDCELL_H
#define GRIDCELL_H
#include <GridCellStockHandler.h>
#include <GridCellCohortHandler.h>
#include <ClimateVariablesCalculator.h>
#include <RevisedTerrestrialPlantModel.h>
/** \file GridCell.h
 * \brief the GridCell header file
 */

class GridCell {
public:
    //----------------------------------------------------------------------------------------------
    //Variables
    //----------------------------------------------------------------------------------------------
    /** \brief The handler for the cohorts in this grid cell */
    GridCellCohortHandler GridCellCohorts;
    /** \brief The handler for the stocks in this grid cell */
    GridCellStockHandler GridCellStocks;
    /** \brief The environmental data for this grid cell */
    map<string, vector<double> > CellEnvironment;
    /** \brief Deltas to track changes in biomasses and abundances of cohorts, stocks and environmental biomass pools during ecological processes */
    map<string, map<string, double>> Deltas;
     /** \brief The latitude of this grid cell */
    float Latitude;
    /** \brief The longitude of this grid cell */
    float Longitude;
    /** \brief An instance of the simple random number generator class */
    std::default_random_engine RandomNumberGenerator;
    /** \brief  Instance of the class to perform general functions*/
    UtilityFunctions Utilities;
    // Optimal prey body size, as a ratio of predator body size
    double OptimalPreyBodySizeRatio;
    //----------------------------------------------------------------------------------------------
    //Methods
    //----------------------------------------------------------------------------------------------
    
    //----------------------------------------------------------------------------------------------
    //default constructor - just to get C++ compile working at the moment
    GridCell() { ;  }
    //----------------------------------------------------------------------------------------------
    /** \brief Constructor for a grid cell; creates cell and reads in environmental data
    @param latitude The latitude of the grid cell 
    @param latIndex The latitudinal index of the grid cell 
    @param longitude The longitude of the grid cell 
    @param lonIndex The longitudinal index of the grid cell 
    @param latCellSize The latitudinal dimension of the grid cell 
    @param lonCellSize The longitudinal dimension of the grid cell 
    @param dataLayers A list of environmental data variables in the model 
    @param missingValue The missing value to be applied to all data in the grid cell 
    @param cohortFunctionalGroups The definitions for cohort functional groups in the model 
    @param stockFunctionalGroups The definitions for stock functional groups in the model 
    @param globalDiagnostics A list of global diagnostic variables for the model grid 
    @param nextCohortID The unique ID number to assign to the next cohort created 
    @param tracking Whether process-tracking is enabled 
    @param specificLocations Whether the model is being run for specific locations */
    GridCell(float latitude, unsigned latIndex, float longitude, unsigned lonIndex, float latCellSize, float lonCellSize,
            map<string, EnviroData>& dataLayers, double missingValue, FunctionalGroupDefinitions& cohortFunctionalGroups,
            FunctionalGroupDefinitions& stockFunctionalGroups, map<string, double>& globalDiagnostics, bool tracking,
            bool specificLocations) {

        // bool to track when environmental data are missing
        bool EnviroMissingValue;

        // Temporary vector for holding initial values of grid cell properties
        vector<double> tempVector;

        // Initialize delta abundance sorted list with appropriate processes
        map<string, double> DeltaAbundance;
        DeltaAbundance["mortality"] = 0.0;

        // Add delta abundance sorted list to deltas sorted list
        Deltas["abundance"] = DeltaAbundance;

        // Initialize delta biomass sorted list with appropriate processes
        map<string, double> DeltaBiomass;
        DeltaBiomass["metabolism"] = 0.0;
        DeltaBiomass["predation"] = 0.0;
        DeltaBiomass["herbivory"] = 0.0;
        DeltaBiomass["reproduction"] = 0.0;

        // Add delta biomass sorted list to deltas sorted list
        Deltas["biomass"] = DeltaBiomass;

        // Initialize delta reproductive biomass vector with appropriate processes
        map<string, double> DeltaReproductiveBiomass;
        DeltaReproductiveBiomass["reproduction"] = 0.0;

        // Add delta reproduction sorted list to deltas sorted list
        Deltas["reproductivebiomass"] = DeltaReproductiveBiomass;

        // Initialize organic pool delta vector with appropriate processes
        map<string, double> DeltaOrganicPool;
        DeltaOrganicPool["herbivory"] = 0.0;
        DeltaOrganicPool["predation"] = 0.0;
        DeltaOrganicPool["mortality"] = 0.0;

        // Add delta organic pool sorted list to deltas sorted list
        Deltas["organicpool"] = DeltaOrganicPool;

        // Initialize respiratory CO2 pool delta vector with appropriate processes
        map<string, double> DeltaRespiratoryCO2Pool;
        DeltaRespiratoryCO2Pool["metabolism"] = 0.0;

        // Add delta respiratory CO2 pool to deltas sorted list
        Deltas["respiratoryCO2pool"] = DeltaRespiratoryCO2Pool;

        // Set the grid cell values of latitude, longitude and missing value as specified
        Latitude = latitude;
        Longitude = longitude;

        //Add the latitude and longitude

        CellEnvironment["Latitude"].push_back(latitude);
        CellEnvironment["Longitude"].push_back(longitude);

        // Add an organic matter pool to the cell environment to track organic biomass not held by animals or plants with an initial value of 0
        CellEnvironment["Organic Pool"].push_back(0.0);

        // Add a repsiratory CO2 pool to the cell environment with an initial value of 0
        CellEnvironment["Respiratory CO2 Pool"].push_back(0.0);

        // Add the grid cell area (in km2) to the cell environment with an initial value of 0
        // Calculate the area of this grid cell
        // Add it to the cell environment
        CellEnvironment["Cell Area"].push_back(Utilities.CalculateGridCellArea(latitude, lonCellSize, latCellSize));

        //Add the latitude and longitude indices
        CellEnvironment["LatIndex"].push_back(latIndex);
        CellEnvironment["LonIndex"].push_back(lonIndex);


        // Add the missing value of data in the grid cell to the cell environment

        CellEnvironment["Missing Value"].push_back(missingValue);

        // Loop through environmental data layers and extract values for this grid cell
        // Also standardise missing values

        // Loop over variables in the list of environmental data
        for (auto Layer : dataLayers) {
            // Initialise the temporary vector of values to be equal to the number of time intervals in the environmental variable
            vector<double> tempVector(Layer.second.NumTimes);
            //Loop over the time intervals in the environmental variable
            for (int hh = 0; hh < Layer.second.NumTimes; hh++) {
                // Add the value of the environmental variable at this time interval to the temporary vector
                tempVector[hh] = Layer.second.GetValue(Latitude, Longitude, (unsigned) hh, EnviroMissingValue, latCellSize, lonCellSize);
                // If the environmental variable is a missing value, then change the value to equal the standard missing value for this cell
                if (EnviroMissingValue)
                    tempVector[hh] = missingValue;
            }
            // Add the values of the environmental variables to the cell environment, with the name of the variable as the key
            CellEnvironment[Layer.first] = tempVector;
        }

        if (CellEnvironment.count("LandSeaMask") != 0) {
            if (CellEnvironment["LandSeaMask"][0] == 0) {
                if (ContainsData(CellEnvironment["OceanTemp"], CellEnvironment["Missing Value"][0])) {
                    //This is a marine cell
                    CellEnvironment["Realm"].push_back(2.0);

                    CellEnvironment["NPP"] = CellEnvironment["OceanNPP"];
                    CellEnvironment["DiurnalTemperatureRange"] = CellEnvironment["OceanDTR"];
                    if (CellEnvironment.count("Temperature") != 0) {
                        if (CellEnvironment.count("SST") != 0) {
                            CellEnvironment["Temperature"] = CellEnvironment["SST"];
                        } else {
                        }
                    } else {
                        CellEnvironment["Temperature"] = CellEnvironment["SST"];
                    }

                } else {
                    //This is a freshwater cell and in this model formulation is characterised as belonging to the terrestrial realm

                    CellEnvironment["Realm"].push_back(1.0);

                    CellEnvironment["NPP"] = CellEnvironment["LandNPP"];
                    CellEnvironment["DiurnalTemperatureRange"] = CellEnvironment["LandDTR"];
                }
            } else {
                //This is a land cell

                CellEnvironment["Realm"].push_back(1.0);

                CellEnvironment["NPP"] = CellEnvironment["LandNPP"];
                CellEnvironment["DiurnalTemperatureRange"] = CellEnvironment["LandDTR"];

            }
        } else {
            cout << "No land sea mask defined - a mask is required to initialise appropriate ecology" << endl;
            exit(1);
        }

        //Calculate and add the standard deviation of monthly temperature as a measure of seasonality
        //Also calculate and add the annual mean temperature for this cell

        vector<double> sdtemp(12);
        vector<double> meantemp(12);

        vector<double> tempTVector = CellEnvironment["Temperature"];

        double Average = 0;
        for (double d : tempTVector)Average += d;
        Average = Average / tempTVector.size();
        meantemp[0] = Average;
        double SumOfSquaresDifferences = 0;
        for (double d : tempTVector)SumOfSquaresDifferences += (d - Average) * (d - Average);
        sdtemp[0] = sqrt(SumOfSquaresDifferences / tempTVector.size());

        CellEnvironment["SDTemperature"] = sdtemp;
        CellEnvironment["AnnualTemperature"] = meantemp;

        //Remove unrequired cell environment layers
        if (CellEnvironment.count("LandNPP") != 0) CellEnvironment.erase("LandNPP");
        if (CellEnvironment.count("LandDTR") != 0) CellEnvironment.erase("LandDTR");
        if (CellEnvironment.count("OceanNPP") != 0) CellEnvironment.erase("OceanNPP");
        if (CellEnvironment.count("OceanDTR") != 0) CellEnvironment.erase("OceanDTR");
        if (CellEnvironment.count("SST") != 0) CellEnvironment.erase("SST");

        // CREATE NPP SEASONALITY LAYER
        CellEnvironment["Seasonality"] = CalculateNPPSeasonality(CellEnvironment["NPP"], CellEnvironment["Missing Value"][0]);

        // Calculate other climate variables from temperature and precipitation
        // Declare an instance of the climate variables calculator
        ClimateVariablesCalculator CVC;

        // Calculate the fraction of the year that experiences frost

        double NDF = CVC.GetNDF(CellEnvironment["FrostDays"], CellEnvironment["Temperature"], CellEnvironment["Missing Value"][0]);
        CellEnvironment["Fraction Year Frost"].push_back(NDF);

        vector<double> frostMonthly(12);
        frostMonthly[0] = min(CellEnvironment["FrostDays"][0] / 31.0, 1.0);
        frostMonthly[1] = min(CellEnvironment["FrostDays"][1] / 28.0, 1.0);
        frostMonthly[2] = min(CellEnvironment["FrostDays"][2] / 31.0, 1.0);
        frostMonthly[3] = min(CellEnvironment["FrostDays"][3] / 30.0, 1.0);
        frostMonthly[4] = min(CellEnvironment["FrostDays"][4] / 31.0, 1.0);
        frostMonthly[5] = min(CellEnvironment["FrostDays"][5] / 30.0, 1.0);
        frostMonthly[6] = min(CellEnvironment["FrostDays"][6] / 31.0, 1.0);
        frostMonthly[7] = min(CellEnvironment["FrostDays"][7] / 31.0, 1.0);
        frostMonthly[8] = min(CellEnvironment["FrostDays"][8] / 30.0, 1.0);
        frostMonthly[9] = min(CellEnvironment["FrostDays"][9] / 31.0, 1.0);
        frostMonthly[10] = min(CellEnvironment["FrostDays"][10] / 30.0, 1.0);
        frostMonthly[11] = min(CellEnvironment["FrostDays"][11] / 31.0, 1.0);

        CellEnvironment["Fraction Month Frost"] = frostMonthly;
        CellEnvironment.erase("FrostDays");

        // Calculate AET and the fractional length of the fire season
        tuple<vector<double>, double, double> TempTuple = CVC.MonthlyActualEvapotranspirationSoilMoisture(CellEnvironment["AWC"][0], CellEnvironment["Precipitation"], CellEnvironment["Temperature"]);
        CellEnvironment["AET"] = get<0>(TempTuple);
        CellEnvironment["Fraction Year Fire"].push_back(get<2> (TempTuple) / 360);

        // Designate a breeding season for this grid cell, where a month is considered to be part of the breeding season if its NPP is at
        // least 80% of the maximum NPP throughout the whole year
        vector<double> BreedingSeason(12);
        for (int i = 0; i < 12; i++) {
            if (CellEnvironment["Seasonality"][i] / *max_element(CellEnvironment["Seasonality"].begin(), CellEnvironment["Seasonality"].end()) > 0.5) {
                BreedingSeason[i] = 1.0;
            } else {
                BreedingSeason[i] = 0.0;
            }
        }
        CellEnvironment["Breeding Season"] = BreedingSeason;

        // Initialise the grid cell cohort and stock handlers
        GridCellCohorts.setSize(cohortFunctionalGroups.GetNumberOfFunctionalGroups());
        GridCellStocks.setSize(stockFunctionalGroups.GetNumberOfFunctionalGroups());

    }
    //----------------------------------------------------------------------------------------------
    /** \brief Converts any missing values to zeroes
    @param data the data vector to convert 
    @param missingValue Missing data value to be converted to zero 
    @return The data vector with any missing data values converted to zero*/
    vector<double> ConvertMissingValuesToZero(vector<double> data, double missingValue) {
    vector<double> TempArray = data;

        for (int ii = 0; ii < TempArray.size(); ii++) {
            TempArray[ii] = (TempArray[ii] != missingValue) ? TempArray[ii] : 0.0;
        }

        return TempArray;
    }
    //----------------------------------------------------------------------------------------------
    /** \brief Checks if any non-missing value data exists in the vector data
    @param data The data vector to be checked 
    @param missingValue The missing value to which the data will be compared 
    @return True if non missing values are found, false if not*/
    bool ContainsData(vector<double> data, double missingValue) {
        bool ContainsData = false;
        for (int ii = 0; ii < data.size(); ii++) {
            if (data[ii] != missingValue) ContainsData = true;
        }
        return ContainsData;
    }
    //----------------------------------------------------------------------------------------------
    /** \brief Checks if any non-missing value data exists in the vector data 
    @param data The data vector to be checked 
    @param missingValue The missing value to which the data will be compared 
    @return True if non missing values are found, false if not*/
    bool ContainsMissingValue(vector<double> data, double missingValue) {
        bool ContainsMV = false;
        for (int ii = 0; ii < data.size(); ii++) {
            if (data[ii] == missingValue) ContainsMV = true;
        }
        return ContainsMV;
    }
    //----------------------------------------------------------------------------------------------
    /** \brief Calculate monthly seasonality values of Net Primary Production - ignores missing values. If there is no NPP data (ie all zero or missing values)
    then assign 1/12 for each month.
    @param NPP Monthly values of NPP 
    @param missingValue Missing data value to which the data will be compared against 
    @return The contribution that each month's NPP makes to annual NPP*/
    vector<double> CalculateNPPSeasonality(vector<double>& NPP, double missingValue) {

        // Check that the NPP data is of monthly temporal resolution
        assert(NPP.size() == 12 && "Error: currently NPP data must be of monthly temporal resolution");

        // Temporary vector to hold seasonality values
        vector<double> NPPSeasonalityValues(12);

        // Loop over months and calculate total annual NPP
        double TotalNPP = 0.0;
        for (int i = 0; i < 12; i++) {
            if (NPP[i] != missingValue && NPP[i] > 0) TotalNPP += NPP[i];
        }
        if (TotalNPP == 0) {
            // Loop over months and calculate seasonality
            // If there is no NPP value then asign a uniform flat seasonality
            for (int i = 0; i < 12; i++) {
                NPPSeasonalityValues[i] = 1.0 / 12.0;
            }

        } else {
            // Some NPP data exists for this grid cell so use that to infer the NPP seasonality
            // Loop over months and calculate seasonality
            for (int i = 0; i < 12; i++) {
                if (NPP[i] != missingValue && NPP[i] > 0) {
                    NPPSeasonalityValues[i] = NPP[i] / TotalNPP;
                } else {
                    NPPSeasonalityValues[i] = 0.0;
                }
            }
        }

        return NPPSeasonalityValues;
    }
    //----------------------------------------------------------------------------------------------
    /** \brief Seed initial stocks and cohorts for this grid cell
    @param cohortFunctionalGroups The functional group definitions for cohorts in the model 
    @param stockFunctionalGroups The functional group definitions for stocks in the model 
    @param globalDiagnostics A list of global diagnostic variables 
    @param nextCohortID The ID number to be assigned to the next produced cohort 
    @param tracking boolean to indicate if cohorts are to be tracked in this model 
    @param totalCellTerrestrialCohorts The total number of cohorts to be seeded in each terrestrial grid cell 
    @param totalCellMarineCohorts The total number of cohorts to be seeded in each marine grid cell 
    @param DrawRandomly Whether the model is set to use random draws 
    @param ZeroAbundance Set this parameter to 'true' if you want to seed the cohorts with zero abundance 
     */
    void SeedGridCellCohortsAndStocks(FunctionalGroupDefinitions& cohortFunctionalGroups, FunctionalGroupDefinitions& stockFunctionalGroups,
            map<string, double>& globalDiagnostics, long long& nextCohortID, bool tracking, double totalCellTerrestrialCohorts,
            double totalCellMarineCohorts, bool DrawRandomly, bool ZeroAbundance) {
        SeedGridCellCohorts(cohortFunctionalGroups, CellEnvironment, globalDiagnostics, nextCohortID, tracking,
                totalCellTerrestrialCohorts, totalCellMarineCohorts, DrawRandomly, ZeroAbundance);
        SeedGridCellStocks(stockFunctionalGroups, CellEnvironment, globalDiagnostics);
    }
    //----------------------------------------------------------------------------------------------
    /** \brief Gets the value in this grid cell of a specified environmental variable at a specified time interval
    @param variableName The name of the environmental layer from which to extract the value 
    @param timeInterval The index of the time interval to return data for (i.e. 0 if it is a yearly variable
    or the month index - 0=Jan, 1=Feb etc. - for monthly variables) 
    @param variableFound Returns whether the variable was found in the cell environment 
    @return The value in this grid cell of a specified environmental variable at a specified time interval*/
    double GetEnviroLayer(string variableName, unsigned timeInterval, bool variableFound) {
        // If the specified variable is in the cell environment then return the requested value, otherwise set variable found boolean
        // to false and return a missing value
        if (CellEnvironment.count(variableName) != 0) {
            variableFound = true;
            return CellEnvironment[variableName][timeInterval];
        } else {
            variableFound = false;
            cout << "Attempt to get environmental layer value failed: " << variableName << " does not exist" << endl;
            return CellEnvironment["Missing Value"][0];
        }
    }
    //----------------------------------------------------------------------------------------------
    /** \brief Sets the value in this grid cell of a specified environmental variable at a specified time interval
    @param variableName The name of the environmental layer to set the value for 
    @param timeInterval The index of the time interval to return data for (i.e. 0 if it is a yearly variable
    or the month index - 0=Jan, 1=Feb etc. - for monthly variables) 
    @param setValue Value to set 
    @return Whether the variable was found in the cell environment*/
    bool SetEnviroLayer(string variableName, unsigned timeInterval, double setValue) {
        // If the specified variable exists in the cell environment then set the specified value and return true; otherwise print an error message and return false
        if (CellEnvironment.count(variableName)) {
            CellEnvironment[variableName][timeInterval] = setValue;
            return true;
        } else {
            cout << "Attempt to set environmental layer value failed: " << variableName << " does not exist" << endl;
            return false;
        }
    }
    //----------------------------------------------------------------------------------------------
    /** \brief Sets the value in this grid cell of a delta of specified type and for a specified ecological process
    @param deltaType The type of delta to set the value for: 'biomass', 'abundance', 'reproductivebiomass', 'organicpool' or 'respiratoryCO2pool 
    @param ecologicalProcess The ecological process to set the value for 
    @param setValue Value to set 
    @return Whether the delta type and ecological process were found within the cell deltas*/
    bool SetDelta(string deltaType, string ecologicalProcess, double setValue) {
        // If the specified ecological and process exist in the cell deltas, then set the value and return true; otherwise, return false
        if (Deltas.count(deltaType) != 0) {
            if (Deltas[deltaType].count(ecologicalProcess) != 0) {
                Deltas[deltaType][ecologicalProcess] = setValue;
                return true;
            } else {
                cout << "Attempt to set delta failed: ecological process " << ecologicalProcess << " does not exist in the list" << endl;
                return false;
            }
        } else {
            cout << "Attempt to set delta failed: delta type " << deltaType << " does not exist in the list" << endl;
            return false;
        }
    }
    //----------------------------------------------------------------------------------------------
    /** \brief  Seed grid cell with cohorts, as specified in the model input files
    @param functionalGroups The functional group definitions for cohorts in the grid cell 
    @param cellEnvironment The environment in the grid cell 
    @param globalDiagnostics A list of global diagnostic variables 
    @param nextCohortID YThe unique ID to assign to the next cohort produced 
    @param tracking boolean to indicate if cohorts are to be tracked in this model 
    @param totalCellTerrestrialCohorts The total number of cohorts to be seeded in each terrestrial grid cell 
    @param totalCellMarineCohorts The total number of cohorts to be seeded in each marine grid cell 
    @param drawRandomly Whether the model is set to use random draws 
    @param ZeroAbundance Set this parameter to 'true' if you want to seed the cohorts with zero abundance */
    void SeedGridCellCohorts(FunctionalGroupDefinitions& functionalGroups, map<string, vector<double>>&
            cellEnvironment, map<string, double>& globalDiagnostics, long long& nextCohortID, bool tracking, double& totalCellTerrestrialCohorts,
            double& totalCellMarineCohorts, bool DrawRandomly, bool ZeroAbundance) {
        // Set the seed for the random number generator from the system time
        unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
        RandomNumberGenerator.seed(seed);
        //
        //            // StreamWriter tempsw = new StreamWriter("C://Temp//adult_juvenile_masses.txt");
        //            // tempsw.WriteLine("adult mass\tjuvenilemass");

        // Define local variables
        double CohortJuvenileMass;
        double CohortAdultMassRatio;
        double CohortAdultMass;
        double ExpectedLnAdultMassRatio;
        vector<int> FunctionalGroupsToUse;
        double NumCohortsThisCell;
        double TotalNewBiomass = 0.0;

        // Get the minimum and maximum possible body masses for organisms in each functional group
        vector<double> MassMinima = functionalGroups.GetBiologicalPropertyAllFunctionalGroups("minimum mass");
        vector<double> MassMaxima = functionalGroups.GetBiologicalPropertyAllFunctionalGroups("maximum mass");

        vector<double> ProportionTimeActive = functionalGroups.GetBiologicalPropertyAllFunctionalGroups("proportion suitable time active");

        //Variable for altering the juvenile to adult mass ratio for marine cells when handling certain functional groups eg baleen whales
        double Scaling = 0.0;

        // Check which realm the cell is in
        if (cellEnvironment["Realm"][0] == 1.0) {

            // Get the indices of all terrestrial functional groups 
            FunctionalGroupsToUse = functionalGroups.GetFunctionalGroupIndex("realm", "terrestrial", true);
            NumCohortsThisCell = totalCellTerrestrialCohorts;
        } else {
            // Get the indices of all marine functional groups

            FunctionalGroupsToUse = functionalGroups.GetFunctionalGroupIndex("realm", "marine", true);
            NumCohortsThisCell = totalCellMarineCohorts;

        }
        assert(cellEnvironment["Realm"][0] > 0.0 && "Missing realm for grid cell");

        if (NumCohortsThisCell > 0);
        {
            //Loop over all functional groups in the model
            for (int FunctionalGroup = 0; FunctionalGroup < functionalGroups.GetNumberOfFunctionalGroups(); FunctionalGroup++) {
                // If it is a functional group that corresponds to the current realm, then seed cohorts
                if (find(FunctionalGroupsToUse.begin(), FunctionalGroupsToUse.end(), FunctionalGroup) != FunctionalGroupsToUse.end()) {

                    // Loop over the initial number of cohorts
                    double NumberOfCohortsInThisFunctionalGroup = 1.0;
                    if (!ZeroAbundance) {
                        NumberOfCohortsInThisFunctionalGroup = functionalGroups.GetBiologicalPropertyOneFunctionalGroup("initial number of gridcellcohorts", FunctionalGroup);
                    }
                    for (int jj = 0; jj < NumberOfCohortsInThisFunctionalGroup; jj++) {
                        // Check whether the model is set to randomly draw the body masses of new cohorts
                        if (DrawRandomly) {
                            // Draw adult mass from a log-normal distribution with mean -6.9 and standard deviation 10.0,
                            // within the bounds of the minimum and maximum body masses for the functional group
                            std::uniform_real_distribution<double> randomNumber(0.0, 1.0);
                            CohortAdultMass = pow(10, (randomNumber(RandomNumberGenerator) * (log10(MassMaxima[FunctionalGroup]) - log10(50 * MassMinima[FunctionalGroup])) + log10(50 * MassMinima[FunctionalGroup])));

                            // Terrestrial and marine organisms have different optimal prey/predator body mass ratios
                            if (cellEnvironment["Realm"][0] == 1.0) {
                                // Optimal prey body size 10%
                                std::normal_distribution<double> randomNumber(0.1, 0.02);
                                OptimalPreyBodySizeRatio = max(0.01, randomNumber(RandomNumberGenerator));
                            } else {
                                if (functionalGroups.GetTraitNames("Diet", FunctionalGroup) == "allspecial") {
                                    // Note that for this group
                                    // it is actually (despite the name) not an optimal prey body size ratio, but an actual body size.
                                    // This is because it is invariant as the predator (filter-feeding baleen whale) grows.
                                    // See also the predation classes.
                                    std::normal_distribution<double> randomNumber(0.0001, 0.1);
                                    OptimalPreyBodySizeRatio = max(0.00001, randomNumber(RandomNumberGenerator));
                                } else {
                                    // Optimal prey body size or marine organisms is 10%
                                    std::normal_distribution<double> randomNumber(0.1, 0.02);
                                    OptimalPreyBodySizeRatio = max(0.01, randomNumber(RandomNumberGenerator));
                                }

                            }

                            // Draw from a log-normal distribution with mean 10.0 and standard deviation 5.0, then add one to obtain 
                            // the ratio of adult to juvenile body mass, and then calculate juvenile mass based on this ratio and within the
                            // bounds of the minimum and maximum body masses for this functional group
                            if (cellEnvironment["Realm"][0] == 1.0) {
                                do {
                                    ExpectedLnAdultMassRatio = 2.24 + 0.13 * log(CohortAdultMass);
                                    std::lognormal_distribution<double> randomNumber(ExpectedLnAdultMassRatio, 0.5);
                                    CohortAdultMassRatio = 1.0 + randomNumber(RandomNumberGenerator);
                                    CohortJuvenileMass = CohortAdultMass * 1.0 / CohortAdultMassRatio;
                                } while (CohortAdultMass <= CohortJuvenileMass || CohortJuvenileMass < MassMinima[FunctionalGroup]);
                            }                                // In the marine realm, have a greater difference between the adult and juvenile body masses, on average
                            else {
                                unsigned Counter = 0;
                                Scaling = 0.2;
                                // Use the scaling to deal with baleen whales not having such a great difference
                                do {

                                    ExpectedLnAdultMassRatio = 2.5 + Scaling * log(CohortAdultMass);
                                    std::lognormal_distribution<double> randomNumber(ExpectedLnAdultMassRatio, 0.5);
                                    CohortAdultMassRatio = 1.0 + 10 * randomNumber(RandomNumberGenerator);
                                    CohortJuvenileMass = CohortAdultMass * 1.0 / CohortAdultMassRatio;
                                    Counter++;
                                    if (Counter > 10) {
                                        Scaling -= 0.01;
                                        Counter = 0;
                                    }
                                } while (CohortAdultMass <= CohortJuvenileMass || CohortJuvenileMass < MassMinima[FunctionalGroup]);
                            }
                        } else {
                            // Use the same seed for the random number generator every time
                            RandomNumberGenerator.seed((unsigned) (jj + 1));

                            // Draw adult mass from a log-normal distribution with mean -6.9 and standard deviation 10.0,
                            // within the bounds of the minimum and maximum body masses for the functional group
                            std::uniform_real_distribution<double> randomNumber(0.0, 1.0);
                            CohortAdultMass = pow(10, (randomNumber(RandomNumberGenerator) * (log10(MassMaxima[FunctionalGroup]) - log10(50 * MassMinima[FunctionalGroup])) + log10(50 * MassMinima[FunctionalGroup])));
                            std::normal_distribution<double> randomNumberN(0.1, 0.02);
                            OptimalPreyBodySizeRatio = max(0.01, randomNumberN(RandomNumberGenerator));

                            // Draw from a log-normal distribution with mean 10.0 and standard deviation 5.0, then add one to obtain 
                            // the ratio of adult to juvenile body mass, and then calculate juvenile mass based on this ratio and within the
                            // bounds of the minimum and maximum body masses for this functional group
                            if (cellEnvironment["Realm"][0] == 1.0) {
                                do {
                                    std::lognormal_distribution<double> randomNumber(ExpectedLnAdultMassRatio, 0.5);
                                    ExpectedLnAdultMassRatio = 2.24 + 0.13 * log(CohortAdultMass);
                                    CohortAdultMassRatio = 1.0 + randomNumber(RandomNumberGenerator);
                                    CohortJuvenileMass = CohortAdultMass * 1.0 / CohortAdultMassRatio;
                                } while (CohortAdultMass <= CohortJuvenileMass || CohortJuvenileMass < MassMinima[FunctionalGroup]);

                            }                                //                                // In the marine realm, have a greater difference between the adult and juvenile body masses, on average
                            else {
                                do {
                                    std::lognormal_distribution<double> randomNumber(ExpectedLnAdultMassRatio, 0.5);
                                    ExpectedLnAdultMassRatio = 2.24 + 0.13 * log(CohortAdultMass);
                                    CohortAdultMassRatio = 1.0 + 10 * randomNumber(RandomNumberGenerator);
                                    ;
                                    CohortJuvenileMass = CohortAdultMass * 1.0 / CohortAdultMassRatio;
                                } while (CohortAdultMass <= CohortJuvenileMass || CohortJuvenileMass < MassMinima[FunctionalGroup]);
                            }
                        }
                        //MB - commented out in original?
                        //double NewBiomass = Math.Pow(0.2, (Math.Log10(CohortAdultMass))) * (1.0E9 * CellEnvironment["Cell Area"][0]) / NumCohortsThisCell;
                        // 3000*(0.6^log(mass)) gives individual cohort biomass density in g ha-1
                        // * 100 to give g km-2
                        // * cell area to give g grid cell
                        //*3300/NumCohortsThisCell scales total initial biomass in the cell to some approximately reasonable mass
                        double NewBiomass = (3300 / NumCohortsThisCell) * 100 * 3000 *
                                pow(0.6, (log10(CohortJuvenileMass))) * (CellEnvironment["Cell Area"][0]);
                        TotalNewBiomass += NewBiomass;
                        double NewAbund = 0.0;
                        if (!ZeroAbundance) {
                            NewAbund = NewBiomass / CohortJuvenileMass;
                        }

                        // Initialise the new cohort with the relevant properties
                        // An instance of Cohort to hold the new cohort
                        Cohort NewCohort(FunctionalGroup, CohortJuvenileMass, CohortAdultMass, CohortJuvenileMass, NewAbund,
                                OptimalPreyBodySizeRatio, 0, ProportionTimeActive[FunctionalGroup], nextCohortID, tracking);

                        // Add the new cohort to the list of grid cell cohorts
                        GridCellCohorts[FunctionalGroup].push_back(NewCohort);


                        // Increment the variable tracking the total number of cohorts in the model
                        globalDiagnostics["NumberOfCohortsInModel"]++;



                    }

                }
            }
        }
    }
    //----------------------------------------------------------------------------------------------
    /** \brief    Seed grid cell with stocks, as specified in the model input files
    @param functionalGroups A reference to the stock functional group handler 
    @param cellEnvironment The environment in the grid cell 
    @param globalDiagnostics A list of global diagnostic variables for the model grid */
    void SeedGridCellStocks(FunctionalGroupDefinitions& functionalGroups, map<string, vector<double>> &
            cellEnvironment, map<string, double>& globalDiagnostics) {
        // Set the seed for the random number generator from the system time
        unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
        RandomNumberGenerator.seed(seed);
        Stock NewStock;

        // Define local variables
        vector<int> FunctionalGroupsToUse;

        // Get the individual body masses for organisms in each stock functional group
        vector<double> IndividualMass = functionalGroups.GetBiologicalPropertyAllFunctionalGroups("individual mass");

        // Check which realm the cell is in
        if (cellEnvironment["Realm"][0] == 1.0 && CellEnvironment["Precipitation"][0] != CellEnvironment["Missing Value"][0] && CellEnvironment["Temperature"][0] != CellEnvironment["Missing Value"][0]) {
            // Get the indices of all terrestrial functional groups 
            FunctionalGroupsToUse = functionalGroups.GetFunctionalGroupIndex("realm", "terrestrial", true);
        } else if (cellEnvironment["Realm"][0] == 2.0 && CellEnvironment["NPP"][0] != CellEnvironment["Missing Value"][0]) {
            // Get the indices of all marine functional groups
            FunctionalGroupsToUse = functionalGroups.GetFunctionalGroupIndex("realm", "marine", true);
        } else {
            // For cells without a realm designation, no functional groups will be used
            //FunctionalGroupsToUse = 0;
        }

        // Loop over all functional groups in the model
        for (int FunctionalGroup = 0; FunctionalGroup < functionalGroups.GetNumberOfFunctionalGroups(); FunctionalGroup++) {
            // If it is a functional group that corresponds to the current realm, then seed the stock
            if (find(FunctionalGroupsToUse.begin(), FunctionalGroupsToUse.end(), FunctionalGroup) != FunctionalGroupsToUse.end()) {
                if (CellEnvironment["Realm"][0] == 1.0) {
                    // An instance of the terrestrial carbon model class
                    RevisedTerrestrialPlantModel PlantModel;

                    // Calculate predicted leaf mass at equilibrium for this stock
                    double LeafMass = PlantModel.CalculateEquilibriumLeafMass(CellEnvironment, functionalGroups.GetTraitNames("leaf strategy", FunctionalGroup) == "deciduous");

                    // Initialise the new stock with the relevant properties
                    Stock NewStock(FunctionalGroup, IndividualMass[FunctionalGroup], LeafMass);

                    // Add the new stock to the list of grid cell stocks
                    GridCellStocks[FunctionalGroup].push_back(NewStock);

                    // Increment the variable tracking the total number of stocks in the model
                    globalDiagnostics["NumberOfStocksInModel"]++;


                } else if (find(FunctionalGroupsToUse.begin(), FunctionalGroupsToUse.end(), FunctionalGroup) != FunctionalGroupsToUse.end()) {
                    // Initialise the new stock with the relevant properties
                    Stock NewStock(FunctionalGroup, IndividualMass[FunctionalGroup], 1e12);

                    // Add the new stock to the list of grid cell stocks
                    GridCellStocks[FunctionalGroup].push_back(NewStock);

                    // Increment the variable tracking the total number of stocks in the model
                    globalDiagnostics["NumberOfStocksInModel"]++;

                } else {
                }

            }

        }

    }
    //----------------------------------------------------------------------------------------------
};
#endif
