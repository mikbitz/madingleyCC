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
    //----------------------------------------------------------------------------------------------
    //Methods
    //----------------------------------------------------------------------------------------------
    
    //----------------------------------------------------------------------------------------------
    //default constructor - just to get grid set up initially
    GridCell() {;}
    //----------------------------------------------------------------------------------------------
    void setCellCoords(float latitude, unsigned latIndex, float longitude, unsigned lonIndex, float latCellSize, float lonCellSize,
             double missingValue){
        // set values for this grid cell
        // Also standardise missing values
    
            // Set the grid cell values of latitude, longitude and missing value as specified
        Latitude = latitude;
        Longitude = longitude;

        //Add the latitude and longitude

        CellEnvironment["Latitude"].push_back(latitude);
        CellEnvironment["Longitude"].push_back(longitude);
        
        // Add the grid cell area (in km2) to the cell environment with an initial value of 0
        // Calculate the area of this grid cell
        // Add it to the cell environment
        CellEnvironment["Cell Area"].push_back(Utilities.CalculateGridCellArea(latitude, lonCellSize, latCellSize));

        //Add the latitude and longitude indices
        CellEnvironment["LatIndex"].push_back(latIndex);
        CellEnvironment["LonIndex"].push_back(lonIndex);
    }    
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
    void SetCellValue(float latitude, unsigned latIndex, float longitude, unsigned lonIndex, float latCellSize, float lonCellSize,
             double missingValue) {


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



        // Add the missing value of data in the grid cell to the cell environment

        CellEnvironment["Missing Value"].push_back(missingValue);
        // Add an organic matter pool to the cell environment to track organic biomass not held by animals or plants with an initial value of 0
        CellEnvironment["Organic Pool"].push_back(0.0);

        // Add a repsiratory CO2 pool to the cell environment with an initial value of 0
        CellEnvironment["Respiratory CO2 Pool"].push_back(0.0);

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
        //GridCellCohorts.setSize(cohortFunctionalGroups.GetNumberOfFunctionalGroups());
        //GridCellStocks.setSize(stockFunctionalGroups.GetNumberOfFunctionalGroups());

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
            // If there is no NPP value then assign a uniform flat seasonality
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
    /** \brief Sets the value in this grid cell of a specified environmental variable at a specified time interval
    @param variableName The name of the environmental layer to set the value for 
    @param timeInterval The index of the time interval to return data for (i.e. 0 if it is a yearly variable
    or the month index - 0=Jan, 1=Feb etc. - for monthly variables) 
    @param setValue Value to set 
    @return Whether the variable was found in the cell environment*/
    bool AddToEnviroLayer(string variableName, unsigned timeInterval, double addValue) {
        // If the specified variable exists in the cell environment then set the specified value and return true; otherwise print an error message and return false
        if (CellEnvironment.count(variableName)) {
            CellEnvironment[variableName][timeInterval] += addValue;
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
    //Apply any function to all cohorts in the cell
    template <typename F>
    void ask(F f) {
    // Loop through functional groups, and perform dispersal according to cohort type and status
        for (int  FG=0; FG < GridCellCohorts.size(); FG++) {
            // Work through the list of cohorts - be sure to work backward to get order of deletion right later
            for (Cohort& c : GridCellCohorts[FG]) {
                f(c);
            }
        }
    }
    //----------------------------------------------------------------------------------------------
    double Realm(){
     return CellEnvironment["Realm"][0];
    }
    //----------------------------------------------------------------------------------------------
    bool isMarine(){
        return (CellEnvironment["Realm"][0]==2.0);
    }
};
#endif
