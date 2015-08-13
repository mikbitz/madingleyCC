#ifndef GRIDCELL_H
#define GRIDCELL_H
#include <GridCellStockHandler.h>
#include <GridCellCohortHandler.h>
/** \file GridCell.h
 * \brief the GridCell header file
 */

//
//namespace Madingley
//{
/** \brief Stores properties of grid cells
<todoD>Remove single valued state-variables and convert model to work with functional groups</todoD>
<todoD>Check the get/set methods and overloads</todoD>
<todoD>Convert GetEnviroLayer to field terminology</todoD> */
class GridCell
    {
    public:
/** \brief The handler for the cohorts in this grid cell */
        GridCellCohortHandler GridCellCohorts;
/** \brief The handler for the stocks in this grid cell */
        GridCellStockHandler GridCellStocks;
/** \brief The environmental data for this grid cell */
        map<string, vector<double> > CellEnvironment;
//        // A sorted list of deltas
/** \brief Deltas to track changes in biomasses and abundances of cohorts, stocks and environmental biomass pools during ecological processes */
          map<string, map<string, double>> Deltas;
/** \brief Get the delta biomasses and abundances for this grid cell */

/** \brief The latitude of this grid cell */
        float Latitude;

/** \brief The longitude of this grid cell */
        float Longitude;
/** \brief
Instance of random number generator to take a time-dependent seed
*/
//        private NonStaticRNG RandomNumberGenerator = new NonStaticRNG();
//        
/** \brief
Instance of the class to perform general functions
*/
       UtilityFunctions Utilities;

       // Optimal prey body size, as a ratio of predator body size
        double OptimalPreyBodySizeRatio;

//default constructor - just to get C++ compile working at the moment
GridCell(){;}
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
           map<string, EnviroData> dataLayers, double missingValue, FunctionalGroupDefinitions cohortFunctionalGroups, 
           FunctionalGroupDefinitions stockFunctionalGroups, map<string, double> globalDiagnostics,bool tracking,
           bool specificLocations)
       {
//
//            // bool to track when environmental data are missing
//            bool EnviroMissingValue;
//
//            // Initialise the utility functions
//            Utilities = new UtilityFunctions();
//
//            // Temporary vector for holding initial values of grid cell properties
//            double[] tempVector;
//
//            // Initialise deltas sorted list
//            _Deltas = new Dictionary<string, Dictionary<string, double>>();
//
//            // Initialize delta abundance sorted list with appropriate processes
//            Dictionary<string, double> DeltaAbundance = new Dictionary<string, double>();
//            DeltaAbundance.Add("mortality", 0.0);
//
//            // Add delta abundance sorted list to deltas sorted list
//            _Deltas.Add("abundance", DeltaAbundance);
//
//            // Initialize delta biomass sorted list with appropriate processes
//            Dictionary<string, double> DeltaBiomass = new Dictionary<string, double>();
//            DeltaBiomass.Add("metabolism", 0.0);
//            DeltaBiomass.Add("predation", 0.0);
//            DeltaBiomass.Add("herbivory", 0.0);
//            DeltaBiomass.Add("reproduction", 0.0);
//
//            // Add delta biomass sorted list to deltas sorted list
//            _Deltas.Add("biomass", DeltaBiomass);
//
//            // Initialize delta reproductive biomass vector with appropriate processes
//            Dictionary<string, double> DeltaReproductiveBiomass = new Dictionary<string, double>();
//            DeltaReproductiveBiomass.Add("reproduction", 0.0);
//
//            // Add delta reproduction sorted list to deltas sorted list
//            _Deltas.Add("reproductivebiomass", DeltaReproductiveBiomass);
//
//            // Initialize organic pool delta vector with appropriate processes
//            Dictionary<string, double> DeltaOrganicPool = new Dictionary<string, double>();
//            DeltaOrganicPool.Add("herbivory", 0.0);
//            DeltaOrganicPool.Add("predation", 0.0);
//            DeltaOrganicPool.Add("mortality", 0.0);
//
//            // Add delta organic pool sorted list to deltas sorted list
//            _Deltas.Add("organicpool", DeltaOrganicPool);
//
//            // Initialize respiratory CO2 pool delta vector with appropriate processes
//            Dictionary<string, double> DeltaRespiratoryCO2Pool = new Dictionary<string, double>();
//            DeltaRespiratoryCO2Pool.Add("metabolism", 0.0);
//
//            // Add delta respiratory CO2 pool to deltas sorted list
//            _Deltas.Add("respiratoryCO2pool", DeltaRespiratoryCO2Pool);
//
//            // Set the grid cell values of latitude, longitude and missing value as specified
//            _Latitude = latitude;
//            _Longitude = longitude;
//
//
//            // Initialise list of environmental data layer values
//            _CellEnvironment = new map<string, double[]>();
//
//            //Add the latitude and longitude
//            tempVector = new double[1];
//            tempVector[0] = latitude;
//            _CellEnvironment.Add("Latitude", tempVector);
//            tempVector = new double[1];
//            tempVector[0] = longitude;
//            _CellEnvironment.Add("Longitude", tempVector);
//
//
//            // Add an organic matter pool to the cell environment to track organic biomass not held by animals or plants with an initial value of 0
//            tempVector = new double[1];
//            tempVector[0] = 0.0;
//            _CellEnvironment.Add("Organic Pool", tempVector);
//
//            // Add a repsiratory CO2 pool to the cell environment with an initial value of 0
//            tempVector = new double[1];
//            tempVector[0] = 0.0;
//            _CellEnvironment.Add("Respiratory CO2 Pool", tempVector);
//
//            // Add the grid cell area (in km2) to the cell environment with an initial value of 0
//            tempVector = new double[1];
//            // Calculate the area of this grid cell
//            tempVector[0] = Utilities.CalculateGridCellArea(latitude, lonCellSize, latCellSize);
//            // Add it to the cell environment
//            _CellEnvironment.Add("Cell Area", tempVector);
//
//            //Add the latitude and longitude indices
//            tempVector = new double[1];
//            tempVector[0] = latIndex;
//            _CellEnvironment.Add("LatIndex", tempVector);
//            tempVector = new double[1];
//            tempVector[0] = lonIndex;
//            _CellEnvironment.Add("LonIndex", tempVector);
//
//
//            // Add the missing value of data in the grid cell to the cell environment
//            tempVector = new double[1];
//            tempVector[0] = missingValue;
//            _CellEnvironment.Add("Missing Value", tempVector);
//
//            // Loop through environmental data layers and extract values for this grid cell
//            // Also standardise missing values
//
//            // Loop over variables in the list of environmental data
//            foreach (string LayerName in dataLayers.Keys)
//            {
//                // Initiliase the temporary vector of values to be equal to the number of time intervals in the environmental variable
//                tempVector = new double[dataLayers[LayerName].NumTimes];
//                // Loop over the time intervals in the environmental variable
//                for (int hh = 0; hh < dataLayers[LayerName].NumTimes; hh++)
//                {
//                    // Add the value of the environmental variable at this time interval to the temporary vector
//                    tempVector[hh] = dataLayers[LayerName].GetValue(_Latitude, _Longitude, (unsigned)hh, out EnviroMissingValue,latCellSize,lonCellSize);
//                    // If the environmental variable is a missing value, then change the value to equal the standard missing value for this cell
//                    if (EnviroMissingValue)
//                        tempVector[hh] = missingValue;
//                }
//                // Add the values of the environmental variables to the cell environment, with the name of the variable as the key
//                _CellEnvironment.Add(LayerName, tempVector);
//            }
//
//            if (specificLocations)
//            {
//                double[, ,] temp;
//                double[, ,] dtr;
//                double[, ,] precip;
//                double[, ,] frost;
//                double[, ,] oceanairt;
//                double[] dtr1 = new double[12];
//                double[] temp1 = new double[12];
//                double[] precip1 = new double[12];
//                double[] frost1 = new double[12];
//                double[] oceanairt1 = new double[12];
//
//                //Declare a dataset to perform the fetch
//                var ds = DataSet.Open("msds:memory2");
//
//                ds.AddAxisCells("longitude", "degrees_east", longitude, longitude + lonCellSize, lonCellSize);
//                ds.AddAxisCells("latitude", "degrees_north", latitude, latitude + latCellSize, latCellSize);
//                
//                ds.AddClimatologyAxisMonthly();
//
//                ds.Fetch(ClimateParameter.FC_TEMPERATURE, "airt", dataSource: EnvironmentalDataSource.ANY); //this call will create 2D variable on dimensions records and months and fill it with a FetchClimate
//                ds.Fetch(ClimateParameter.FC_LAND_DIURNAL_TEMPERATURE_RANGE, "landDTR", dataSource: EnvironmentalDataSource.ANY); //this call will create 2D variable on dimensions records and months and fill it with a FetchClimate
//                ds.Fetch(ClimateParameter.FC_PRECIPITATION, "precip", dataSource: EnvironmentalDataSource.ANY); //this call will create 2D variable on dimensions records and months and fill it with a FetchClimate
//                ds.Fetch(ClimateParameter.FC_LAND_FROST_DAY_FREQUENCY, "frost", dataSource: EnvironmentalDataSource.ANY);
//                ds.Fetch(ClimateParameter.FC_OCEAN_AIR_TEMPERATURE, "oceanairt", dataSource: EnvironmentalDataSource.ANY);
//
//                dtr = (double[, ,])ds.Variables["landDTR"].GetData();
//                temp = (double[,,])ds.Variables["airt"].GetData();
//                precip = (double[,,])ds.Variables["precip"].GetData();
//                frost = (double[, ,])ds.Variables["frost"].GetData();
//                oceanairt = (double[, ,])ds.Variables["oceanairt"].GetData();
//
//                int m = 0;
//                for (m = 0; m < 12; m++)
//                {
//                    dtr1[m] = dtr[m, 0, 0];
//                    temp1[m] = temp[m, 0, 0];
//                    precip1[m] = precip[m, 0, 0];
//                    frost1[m] = frost[m, 0, 0];
//                    oceanairt1[m] = oceanairt[m, 0, 0];
//                }
//
//
//
//                _CellEnvironment.Add("LandDTR", dtr1);
//                _CellEnvironment.Add("Temperature", temp1);
//                _CellEnvironment.Add("Precipitation", precip1);
//                _CellEnvironment.Add("FrostDays", frost1);
//                _CellEnvironment.Add("OceanTemp", oceanairt1);
//
//
//
//            }
//
//            if (_CellEnvironment.ContainsKey("LandSeaMask"))
//            {
//                if (_CellEnvironment["LandSeaMask"][0].CompareTo(0.0) == 0)
//                {
//                    if (ContainsData(_CellEnvironment["OceanTemp"], _CellEnvironment["Missing Value"][0]))
//                    {
//                        //This is a marine cell
//                        tempVector = new double[1];
//                        tempVector[0] = 2.0;
//                        _CellEnvironment.Add("Realm", tempVector);
//
//                        _CellEnvironment.Add("NPP", _CellEnvironment["OceanNPP"]);
//                        _CellEnvironment.Add("DiurnalTemperatureRange", _CellEnvironment["OceanDTR"]);
//                        if (_CellEnvironment.ContainsKey("Temperature"))
//                        {
//                            if(_CellEnvironment.ContainsKey("SST"))
//                            {
//                                _CellEnvironment["Temperature"] = _CellEnvironment["SST"];
//                            }
//                            else
//                            {
//                            }
//                        }
//                        else
//                        {
//                            _CellEnvironment.Add("Temperature", _CellEnvironment["SST"]);
//                        }
//
//                    }
//                    else
//                    {
//                        //This is a freshwater cell and in this model formulation is characterised as belonging to the terrestrial realm
//                        tempVector = new double[1];
//                        tempVector[0] = 1.0;
//                        _CellEnvironment.Add("Realm", tempVector);
//
//                        _CellEnvironment.Add("NPP", _CellEnvironment["LandNPP"]);
//                        _CellEnvironment.Add("DiurnalTemperatureRange", _CellEnvironment["LandDTR"]);
//                    }
//                }
//                else
//                {
//                    //This is a land cell
//                    tempVector = new double[1];
//                    tempVector[0] = 1.0;
//                    _CellEnvironment.Add("Realm", tempVector);
//
//                    _CellEnvironment.Add("NPP", _CellEnvironment["LandNPP"]);
//                    _CellEnvironment.Add("DiurnalTemperatureRange", _CellEnvironment["LandDTR"]);
//
//                }
//            }
//            else
//            {
//                Debug.Fail("No land sea mask defined - a mask is required to initialise appropriate ecology");
//            }
//
//            //Calculate and add the standard deviation of monthly temperature as a measure of seasonality
//            //Also calculate and add the annual mean temperature for this cell
//            tempVector = new double[12];
//            double[] sdtemp = new double[12];
//            double[] meantemp = new double[12];
//
//            tempVector = _CellEnvironment["Temperature"];
//
//            double Average = tempVector.Average();
//            meantemp[0] = Average;
//            double SumOfSquaresDifferences = tempVector.Select(val => (val - Average) * (val - Average)).Sum();
//            sdtemp[0] = Math.Sqrt(SumOfSquaresDifferences / tempVector.Length);
//
//            _CellEnvironment.Add("SDTemperature", sdtemp);
//            _CellEnvironment.Add("AnnualTemperature", meantemp);
//
//            //Remove unrequired cell environment layers
//            if (_CellEnvironment.ContainsKey("LandNPP")) _CellEnvironment.Remove("LandNPP");
//            if (_CellEnvironment.ContainsKey("LandDTR")) _CellEnvironment.Remove("LandDTR");
//            if (_CellEnvironment.ContainsKey("OceanNPP")) _CellEnvironment.Remove("OceanNPP");
//            if (_CellEnvironment.ContainsKey("OceanDTR")) _CellEnvironment.Remove("OceanDTR");
//            if (_CellEnvironment.ContainsKey("SST")) _CellEnvironment.Remove("SST");
//
//            // CREATE NPP SEASONALITY LAYER
//            _CellEnvironment.Add("Seasonality", CalculateNPPSeasonality(_CellEnvironment["NPP"], _CellEnvironment["Missing Value"][0]));
//
//            // Calculate other climate variables from temperature and precipitation
//            // Declare an instance of the climate variables calculator
//            ClimateVariablesCalculator CVC = new ClimateVariablesCalculator();
//
//            // Calculate the fraction of the year that experiences frost
//            double[] NDF = new double[1];
//            NDF[0] = CVC.GetNDF(_CellEnvironment["FrostDays"], _CellEnvironment["Temperature"],_CellEnvironment["Missing Value"][0]);
//            _CellEnvironment.Add("Fraction Year Frost", NDF);
//
//            double[] frostMonthly = new double[12];
//            frostMonthly[0] = Math.Min(_CellEnvironment["FrostDays"][0] / 31.0, 1.0);
//            frostMonthly[1] = Math.Min(_CellEnvironment["FrostDays"][1] / 28.0, 1.0);
//            frostMonthly[2] = Math.Min(_CellEnvironment["FrostDays"][2] / 31.0, 1.0);
//            frostMonthly[3] = Math.Min(_CellEnvironment["FrostDays"][3] / 30.0, 1.0);
//            frostMonthly[4] = Math.Min(_CellEnvironment["FrostDays"][4] / 31.0, 1.0);
//            frostMonthly[5] = Math.Min(_CellEnvironment["FrostDays"][5] / 30.0, 1.0);
//            frostMonthly[6] = Math.Min(_CellEnvironment["FrostDays"][6] / 31.0, 1.0);
//            frostMonthly[7] = Math.Min(_CellEnvironment["FrostDays"][7] / 31.0, 1.0);
//            frostMonthly[8] = Math.Min(_CellEnvironment["FrostDays"][8] / 30.0, 1.0);
//            frostMonthly[9] = Math.Min(_CellEnvironment["FrostDays"][9] / 31.0, 1.0);
//            frostMonthly[10] = Math.Min(_CellEnvironment["FrostDays"][10] / 30.0, 1.0);
//            frostMonthly[11] = Math.Min(_CellEnvironment["FrostDays"][11] / 31.0, 1.0);
//
//            _CellEnvironment.Add("Fraction Month Frost", frostMonthly);
//            _CellEnvironment.Remove("FrostDays");
//
//            // Calculate AET and the fractional length of the fire season
//            Tuple<double[], double, double> TempTuple = new Tuple<double[], double, double>(new double[12], new double(), new double());
//            TempTuple = CVC.MonthlyActualEvapotranspirationSoilMoisture(_CellEnvironment["AWC"][0], _CellEnvironment["Precipitation"], _CellEnvironment["Temperature"]);
//            _CellEnvironment.Add("AET", TempTuple.Item1);
//            _CellEnvironment.Add("Fraction Year Fire", new double[1] { TempTuple.Item3 / 360 });
//
//            // Designate a breeding season for this grid cell, where a month is considered to be part of the breeding season if its NPP is at
//            // least 80% of the maximum NPP throughout the whole year
//            double[] BreedingSeason = new double[12];
//            for (int i = 0; i < 12; i++)
//            {
//                if ((_CellEnvironment["Seasonality"][i] / _CellEnvironment["Seasonality"].Max()) > 0.5)
//                {
//                    BreedingSeason[i] = 1.0;
//                }
//                else
//                {
//                    BreedingSeason[i] = 0.0;
//                }
//            }
//            _CellEnvironment.Add("Breeding Season", BreedingSeason);
//
//
//            // Initialise the grid cell cohort and stock handlers
//            _GridCellCohorts = new GridCellCohortHandler(cohortFunctionalGroups.GetNumberOfFunctionalGroups());
//            _GridCellStocks = new GridCellStockHandler(stockFunctionalGroups.GetNumberOfFunctionalGroups());

       }

/** \brief Converts any missing values to zeroes

@param data the data vector to convert 
@param missingValue Missing data value to be converted to zero 
@return The data vector with any missing data values converted to zero*/
        vector<double> ConvertMissingValuesToZero(vector<double> data, double missingValue)
       {
           vector<double> TempArray = data;

           for (int ii = 0; ii < TempArray.size(); ii++)
           {
               TempArray[ii] = (TempArray[ii]!=missingValue) ? TempArray[ii] : 0.0;
           }

           return TempArray;
       }

/** \brief Checks if any non-missing value data exists in the vector data

@param data The data vector to be checked 
@param missingValue The missing value to which the data will be compared 
@return True if non missing values are found, false if not*/
        bool ContainsData(vector<double> data, double missingValue)
       {
           bool ContainsData = false;
           for (int ii = 0; ii < data.size(); ii++)
           {
               if (data[ii]!=missingValue) ContainsData = true;
           }
           return ContainsData;
       }



/** \brief Checks if any non-missing value data exists in the vector data 
@param data The data vector to be checked 
@param missingValue The missing value to which the data will be compared 
@return True if non missing values are found, false if not*/
        bool ContainsMissingValue(vector<double>  data, double missingValue)
       {
           bool ContainsMV = false;
           for (int ii = 0; ii < data.size(); ii++)
           {
               if (data[ii]==missingValue) ContainsMV = true;
           }
           return ContainsMV;
       }

/** \brief Calculate monthly seasonality values of Net Primary Production - ignores missing values. If there is no NPP data (ie all zero or missing values)
then assign 1/12 for each month.

@param NPP Monthly values of NPP 
@param missingValue Missing data value to which the data will be compared against 
@return The contribution that each month's NPP makes to annual NPP*/
        vector<double> CalculateNPPSeasonality(vector<double> NPP, double missingValue)
       {

//            // Check that the NPP data is of monthly temporal resolution
//            Debug.Assert(NPP.Length == 12, "Error: currently NPP data must be of monthly temporal resolution");
//
//            // Temporary vector to hold seasonality values
//            double[] NPPSeasonalityValues = new double[12];
//
//            // Loop over months and calculate total annual NPP
//            double TotalNPP = 0.0;
//            for (int i = 0; i < 12; i++)
//            {
//                if (NPP[i].CompareTo(missingValue) != 0 && NPP[i].CompareTo(0.0) > 0) TotalNPP += NPP[i];
//            }
//            if (TotalNPP.CompareTo(0.0) == 0)
//            {
//                // Loop over months and calculate seasonality
//                // If there is no NPP value then asign a uniform flat seasonality
//                for (int i = 0; i < 12; i++)
//                {
//                    NPPSeasonalityValues[i] = 1.0/12.0;
//                }
//
//            }
//            else
//            {
//                // Some NPP data exists for this grid cell so use that to infer the NPP seasonality
//                // Loop over months and calculate seasonality
//                for (int i = 0; i < 12; i++)
//                {
//                    if (NPP[i].CompareTo(missingValue) != 0 && NPP[i].CompareTo(0.0) > 0)
//                    {
//                        NPPSeasonalityValues[i] = NPP[i] / TotalNPP;
//                    }
//                    else
//                    {
//                        NPPSeasonalityValues[i] = 0.0;
//                    }
//                }
//            }
//
//            return NPPSeasonalityValues;
       }

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
        void SeedGridCellCohortsAndStocks(FunctionalGroupDefinitions cohortFunctionalGroups, FunctionalGroupDefinitions stockFunctionalGroups,
           map<string, double> globalDiagnostics, long long& nextCohortID, bool tracking, double totalCellTerrestrialCohorts, 
           double totalCellMarineCohorts, bool DrawRandomly, bool ZeroAbundance)
       {
//            SeedGridCellCohorts(ref cohortFunctionalGroups, ref _CellEnvironment, globalDiagnostics, ref nextCohortID, tracking, 
//                totalCellTerrestrialCohorts, totalCellMarineCohorts, DrawRandomly, ZeroAbundance);
//            SeedGridCellStocks(ref stockFunctionalGroups, ref _CellEnvironment, globalDiagnostics);
       }

/** \brief Gets the value in this grid cell of a specified environmental variable at a specified time interval

@param variableName The name of the environmental layer from which to extract the value 
@param timeInterval The index of the time interval to return data for (i.e. 0 if it is a yearly variable
or the month index - 0=Jan, 1=Feb etc. - for monthly variables) 
@param variableFound Returns whether the variable was found in the cell environment 
@return The value in this grid cell of a specified environmental variable at a specified time interval*/
        double GetEnviroLayer(string variableName, unsigned timeInterval,  bool variableFound)
       {
//            // If the specified variable is in the cell environment then return the requested value, otherwise set variable found boolean
//            // to false and return a missing value
//            if (_CellEnvironment.ContainsKey(variableName))
//            {
//                variableFound = true;
//                return _CellEnvironment[variableName][timeInterval];
//            }
//            else
//            {
//                variableFound = false;
//                Console.WriteLine("Attempt to get environmental layer value failed: {0} does not exist", variableName);
//                return _CellEnvironment["Missing Value"][0];
//            }
       }

/** \brief Sets the value in this grid cell of a specified environmental variable at a specified time interval

@param variableName The name of the environmental layer to set the value for 
@param timeInterval The index of the time interval to return data for (i.e. 0 if it is a yearly variable
or the month index - 0=Jan, 1=Feb etc. - for monthly variables) 
@param setValue Value to set 
@return Whether the variable was found in the cell environment*/
        bool SetEnviroLayer(string variableName, unsigned timeInterval, double setValue)
       {
//            // If the specified variable exists in the cell environment then set the specified value and return true; otherwise print an error message and return false
//            if (_CellEnvironment.ContainsKey(variableName))
//            {
//                _CellEnvironment[variableName][timeInterval] = setValue;
//                return true;
//            }
//            else
//            {
//                Console.WriteLine("Attempt to set environmental layer value failed: {0} does not exist", variableName);
//                return false;
//            }
       }

/** \brief Sets the value in this grid cell of a delta of specified type and for a specified ecological process

@param deltaType The type of delta to set the value for: 'biomass', 'abundance', 'reproductivebiomass', 'organicpool' or 'respiratoryCO2pool 
@param ecologicalProcess The ecological process to set the value for 
@param setValue Value to set 
@return Whether the delta type and ecological process were found within the cell deltas*/
        bool SetDelta(string deltaType, string ecologicalProcess, double setValue)
       {
//            // If the specified ecological and process exist in the cell deltas, then set the value and return true; otherwise, return false
//            if (_Deltas.ContainsKey(deltaType))
//            {
//                if (_Deltas[deltaType].ContainsKey(ecologicalProcess))
//                {
//                    _Deltas[deltaType][ecologicalProcess] = setValue;
//                    return true;
//                }
//                else
//                {
//                    Console.WriteLine("Attempt to set delta failed: ecological process '{0}' does not exist in the list", ecologicalProcess);
//                    return false;
//                }
//            }
//            else
//            {
//                Console.WriteLine("Attempt to set delta failed: delta type '{0}' does not exist in the list", deltaType);
//                return false;
//            }
       }

/** \brief
Seed grid cell with cohorts, as specified in the model input files

@param functionalGroups The functional group definitions for cohorts in the grid cell 
@param cellEnvironment The environment in the grid cell 
@param globalDiagnostics A list of global diagnostic variables 
@param nextCohortID YThe unique ID to assign to the next cohort produced 
@param tracking boolean to indicate if cohorts are to be tracked in this model 
@param totalCellTerrestrialCohorts The total number of cohorts to be seeded in each terrestrial grid cell 
@param totalCellMarineCohorts The total number of cohorts to be seeded in each marine grid cell 
@param drawRandomly Whether the model is set to use random draws 
@param ZeroAbundance Set this parameter to 'true' if you want to seed the cohorts with zero abundance */
//        private void SeedGridCellCohorts(ref FunctionalGroupDefinitions functionalGroups, ref map<string, double[]>
//            cellEnvironment, map<string, double> globalDiagnostics, ref Int64 nextCohortID, bool tracking, double totalCellTerrestrialCohorts, 
//            double totalCellMarineCohorts, bool DrawRandomly, bool ZeroAbundance)
//        {
//            // Set the seed for the random number generator from the system time
//            RandomNumberGenerator.SetSeedFromSystemTime();
//
//            // StreamWriter tempsw = new StreamWriter("C://Temp//adult_juvenile_masses.txt");
//            // tempsw.WriteLine("adult mass\tjuvenilemass");
//
//            // Define local variables
//            double CohortJuvenileMass;
//            double CohortAdultMassRatio;
//            double CohortAdultMass;
//            double ExpectedLnAdultMassRatio;
//            int[] FunctionalGroupsToUse;
//            double NumCohortsThisCell;
//            double TotalNewBiomass =0.0;
//
//            // Get the minimum and maximum possible body masses for organisms in each functional group
//            double[] MassMinima = functionalGroups.GetBiologicalPropertyAllFunctionalGroups("minimum mass");
//            double[] MassMaxima = functionalGroups.GetBiologicalPropertyAllFunctionalGroups("maximum mass");
//
//            double[] ProportionTimeActive = functionalGroups.GetBiologicalPropertyAllFunctionalGroups("proportion suitable time active");
//
//            //Variable for altering the juvenile to adult mass ratio for marine cells when handling certain functional groups eg baleen whales
//            double Scaling = 0.0;
//            
//            // Check which realm the cell is in
//            if (cellEnvironment["Realm"][0] == 1.0)
//            {
//                // Get the indices of all terrestrial functional groups 
//                FunctionalGroupsToUse = functionalGroups.GetFunctionalGroupIndex("realm", "terrestrial", true);
//                NumCohortsThisCell = totalCellTerrestrialCohorts;
//            }
//            else
//            {
//                // Get the indices of all marine functional groups
//                FunctionalGroupsToUse = functionalGroups.GetFunctionalGroupIndex("realm", "marine", true);
//                NumCohortsThisCell = totalCellMarineCohorts;
//            }
//            Debug.Assert(cellEnvironment["Realm"][0] > 0.0, "Missing realm for grid cell");
//
//            if (NumCohortsThisCell > 0)
//            {
//                // Loop over all functional groups in the model
//                for (int FunctionalGroup = 0; FunctionalGroup < functionalGroups.GetNumberOfFunctionalGroups(); FunctionalGroup++)
//                {
//
//                    // Create a new list to hold the cohorts in the grid cell
//                    _GridCellCohorts[FunctionalGroup] = new List<Cohort>();
//
//                    // If it is a functional group that corresponds to the current realm, then seed cohorts
//                    if (FunctionalGroupsToUse.Contains(FunctionalGroup))
//                    {
//                        // Loop over the initial number of cohorts
//                        double NumberOfCohortsInThisFunctionalGroup = 1.0;
//                        if (!ZeroAbundance)
//                        {
//                            NumberOfCohortsInThisFunctionalGroup = functionalGroups.GetBiologicalPropertyOneFunctionalGroup("initial number of gridcellcohorts", FunctionalGroup);
//                        }
//                        for (int jj = 0; jj < NumberOfCohortsInThisFunctionalGroup; jj++)
//                        {
//                            // Check whether the model is set to randomly draw the body masses of new cohorts
//                            if (DrawRandomly)
//                            {
//                                // Draw adult mass from a log-normal distribution with mean -6.9 and standard deviation 10.0,
//                                // within the bounds of the minimum and maximum body masses for the functional group
//                                CohortAdultMass = Math.Pow(10, (RandomNumberGenerator.GetUniform() * (Math.Log10(MassMaxima[FunctionalGroup]) - Math.Log10(50 * MassMinima[FunctionalGroup])) + Math.Log10(50 * MassMinima[FunctionalGroup])));
//
//                                // Terrestrial and marine organisms have different optimal prey/predator body mass ratios
//                                if (cellEnvironment["Realm"][0] == 1.0)
//                                    // Optimal prey body size 10%
//                                    OptimalPreyBodySizeRatio = Math.Max(0.01, RandomNumberGenerator.GetNormal(0.1, 0.02));
//                                else
//                                {
//                                    if (functionalGroups.GetTraitNames("Diet", FunctionalGroup) == "allspecial")
//                                    {
//                                        // Note that for this group
//                                        // it is actually (despite the name) not an optimal prey body size ratio, but an actual body size.
//                                        // This is because it is invariant as the predator (filter-feeding baleen whale) grows.
//                                        // See also the predation classes.
//                                        OptimalPreyBodySizeRatio = Math.Max(0.00001, RandomNumberGenerator.GetNormal(0.0001, 0.1));
//                                    }
//                                    else
//                                    {
//                                        // Optimal prey body size or marine organisms is 10%
//                                        OptimalPreyBodySizeRatio = Math.Max(0.01, RandomNumberGenerator.GetNormal(0.1, 0.02));
//                                    }
//
//                                }
//
//
//                                // Draw from a log-normal distribution with mean 10.0 and standard deviation 5.0, then add one to obtain 
//                                // the ratio of adult to juvenile body mass, and then calculate juvenile mass based on this ratio and within the
//                                // bounds of the minimum and maximum body masses for this functional group
//                                if (cellEnvironment["Realm"][0] == 1.0)
//                                {
//                                    do
//                                    {
//                                        ExpectedLnAdultMassRatio = 2.24 + 0.13 * Math.Log(CohortAdultMass);
//                                        CohortAdultMassRatio = 1.0 + RandomNumberGenerator.GetLogNormal(ExpectedLnAdultMassRatio, 0.5);
//                                        CohortJuvenileMass = CohortAdultMass * 1.0 / CohortAdultMassRatio;
//                                    } while (CohortAdultMass <= CohortJuvenileMass || CohortJuvenileMass < MassMinima[FunctionalGroup]);
//                                }
//                                // In the marine realm, have a greater difference between the adult and juvenile body masses, on average
//                                else
//                                {
//                                    unsigned Counter = 0;
//                                    Scaling = 0.2;
//                                    // Use the scaling to deal with baleen whales not having such a great difference
//                                    do
//                                    {
//
//                                        ExpectedLnAdultMassRatio = 2.5 + Scaling * Math.Log(CohortAdultMass);
//                                        CohortAdultMassRatio = 1.0 + 10 * RandomNumberGenerator.GetLogNormal(ExpectedLnAdultMassRatio, 0.5);
//                                        CohortJuvenileMass = CohortAdultMass * 1.0 / CohortAdultMassRatio;
//                                        Counter++;
//                                        if (Counter > 10)
//                                        {
//                                            Scaling -= 0.01;
//                                            Counter = 0;
//                                        }
//                                    } while (CohortAdultMass <= CohortJuvenileMass || CohortJuvenileMass < MassMinima[FunctionalGroup]);
//                                }
//                            }
//                            else
//                            {
//                                // Use the same seed for the random number generator every time
//                                RandomNumberGenerator.SetSeed((unsigned)(jj + 1));
//
//                                // Draw adult mass from a log-normal distribution with mean -6.9 and standard deviation 10.0,
//                                // within the bounds of the minimum and maximum body masses for the functional group
//                                CohortAdultMass = Math.Pow(10, (RandomNumberGenerator.GetUniform() * (Math.Log10(MassMaxima[FunctionalGroup]) - Math.Log10(50 * MassMinima[FunctionalGroup])) + Math.Log10(50 * MassMinima[FunctionalGroup])));
//                                
//                                OptimalPreyBodySizeRatio = Math.Max(0.01, RandomNumberGenerator.GetNormal(0.1, 0.02));
//                                
//                                // Draw from a log-normal distribution with mean 10.0 and standard deviation 5.0, then add one to obtain 
//                                // the ratio of adult to juvenile body mass, and then calculate juvenile mass based on this ratio and within the
//                                // bounds of the minimum and maximum body masses for this functional group
//                                if (cellEnvironment["Realm"][0] == 1.0)
//                                {
//                                    do
//                                    {
//                                        ExpectedLnAdultMassRatio = 2.24 + 0.13 * Math.Log(CohortAdultMass);
//                                        CohortAdultMassRatio = 1.0 + RandomNumberGenerator.GetLogNormal(ExpectedLnAdultMassRatio, 0.5);
//                                        CohortJuvenileMass = CohortAdultMass * 1.0 / CohortAdultMassRatio;
//                                    } while (CohortAdultMass <= CohortJuvenileMass || CohortJuvenileMass < MassMinima[FunctionalGroup]);
//
//                                }
//                                // In the marine realm, have a greater difference between the adult and juvenile body masses, on average
//                                else
//                                {
//                                    do
//                                    {
//                                        ExpectedLnAdultMassRatio = 2.24 + 0.13 * Math.Log(CohortAdultMass);
//                                        CohortAdultMassRatio = 1.0 + 10 * RandomNumberGenerator.GetLogNormal(ExpectedLnAdultMassRatio, 0.5);
//                                        CohortJuvenileMass = CohortAdultMass * 1.0 / CohortAdultMassRatio;
//                                    } while (CohortAdultMass <= CohortJuvenileMass || CohortJuvenileMass < MassMinima[FunctionalGroup]);
//                                }
//                            }
//
//                            // An instance of Cohort to hold the new cohort
//                            Cohort NewCohort;
//                            //double NewBiomass = Math.Pow(0.2, (Math.Log10(CohortAdultMass))) * (1.0E9 * _CellEnvironment["Cell Area"][0]) / NumCohortsThisCell;
//                            // 3000*(0.6^log(mass)) gives individual cohort biomass density in g ha-1
//                            // * 100 to give g km-2
//                            // * cell area to give g grid cell
//                            //*3300/NumCohortsThisCell scales total initial biomass in the cell to some approximately reasonable mass
//                            double NewBiomass = (3300 / NumCohortsThisCell) * 100 * 3000 * 
//                                Math.Pow(0.6, (Math.Log10(CohortJuvenileMass))) * (_CellEnvironment["Cell Area"][0]);
//                            TotalNewBiomass += NewBiomass;
//                            double NewAbund = 0.0;
//                            if (!ZeroAbundance)
//                            {
//                                NewAbund = NewBiomass / CohortJuvenileMass;
//                            }
//
//                            /*
//                            // TEMPORARILY MARINE ONLY
//                            if (cellEnvironment["Realm"][0] == 1)
//                            {
//                                NewAbund = 0.0;
//                            }
//                            */
//
//                            // Initialise the new cohort with the relevant properties
//                            NewCohort = new Cohort((byte)FunctionalGroup, CohortJuvenileMass, CohortAdultMass, CohortJuvenileMass, NewAbund,
//                            OptimalPreyBodySizeRatio, (ushort)0, ProportionTimeActive[FunctionalGroup], ref nextCohortID, tracking);
//
//                            // Add the new cohort to the list of grid cell cohorts
//                            _GridCellCohorts[FunctionalGroup].Add(NewCohort);
//
//
//
//                            /// TEMPORARY
//                            /*
//                            // Check whether the model is set to randomly draw the body masses of new cohorts
//                            if ((Longitude % 4 == 0) && (Latitude % 4 == 0))
//                            {
//                                if (DrawRandomly)
//                                {
//                                    CohortAdultMass = 100000;
//                                    CohortJuvenileMass = 100000;
//                                }
//                                else
//                                {
//                                    CohortAdultMass = 100000;
//                                    CohortJuvenileMass = 100000;
//                                }
//
//                                // An instance of Cohort to hold the new cohort
//                                Cohort NewCohort;
//                                double NewBiomass = (1.0E7 * _CellEnvironment["Cell Area"][0]) / NumCohortsThisCell;
//                                double NewAbund = 0.0;
//                                NewAbund = 3000;
//
//                                // Initialise the new cohort with the relevant properties
//                                NewCohort = new Cohort((byte)FunctionalGroup, CohortJuvenileMass, CohortAdultMass, CohortJuvenileMass, NewAbund,
//                                    (ushort)0, ref nextCohortID, tracking);
//
//                                // Add the new cohort to the list of grid cell cohorts
//                                _GridCellCohorts[FunctionalGroup].Add(NewCohort);
//                            }
//                            */
//
//                            // Incrememt the variable tracking the total number of cohorts in the model
//                            globalDiagnostics["NumberOfCohortsInModel"]++;
//
//
//
//                        }
//
//                    }
//                }
//
//            }
//            else
//            {
//                                // Loop over all functional groups in the model
//                for (int FunctionalGroup = 0; FunctionalGroup < functionalGroups.GetNumberOfFunctionalGroups(); FunctionalGroup++)
//                {
//
//                    // Create a new list to hold the cohorts in the grid cell
//                    _GridCellCohorts[FunctionalGroup] = new List<Cohort>();
//                }       
//            }
//
//            // tempsw.Dispose();
//        }
//
/** \brief
Seed grid cell with stocks, as specified in the model input files

@param functionalGroups A reference to the stock functional group handler 
@param cellEnvironment The environment in the grid cell 
@param globalDiagnostics A list of global diagnostic variables for the model grid */
//        private void SeedGridCellStocks(ref FunctionalGroupDefinitions functionalGroups, ref map<string, double[]> 
//            cellEnvironment, map<string, double> globalDiagnostics)
//        {
//            // Set the seed for the random number generator from the system time
//            RandomNumberGenerator.SetSeedFromSystemTime();
//
//            Stock NewStock;
//
//            // Define local variables
//            int[] FunctionalGroupsToUse;
//
//            // Get the individual body masses for organisms in each stock functional group
//            double[] IndividualMass = functionalGroups.GetBiologicalPropertyAllFunctionalGroups("individual mass");
//
//            // Check which realm the cell is in
//            if (cellEnvironment["Realm"][0] == 1.0 && _CellEnvironment["Precipitation"][0] != _CellEnvironment["Missing Value"][0] && _CellEnvironment["Temperature"][0] != _CellEnvironment["Missing Value"][0])
//            {
//                // Get the indices of all terrestrial functional groups 
//                FunctionalGroupsToUse = functionalGroups.GetFunctionalGroupIndex("realm", "terrestrial", true);
//            }
//            else if (cellEnvironment["Realm"][0] == 2.0 && _CellEnvironment["NPP"][0] != _CellEnvironment["Missing Value"][0])
//            {
//                // Get the indices of all marine functional groups
//                FunctionalGroupsToUse = functionalGroups.GetFunctionalGroupIndex("realm", "marine", true);
//            }
//            else
//            {
//                // For cells without a realm designation, no functional groups will be used
//                FunctionalGroupsToUse = new int[0];
//            }
//
//            // Loop over all functional groups in the model
//            for (int FunctionalGroup = 0; FunctionalGroup < functionalGroups.GetNumberOfFunctionalGroups(); FunctionalGroup++)
//            {
//                // Create a new list to hold the stocks in the grid cell
//                _GridCellStocks[FunctionalGroup] = new List<Stock>();
//
//                // If it is a functional group that corresponds to the current realm, then seed the stock
//                if (FunctionalGroupsToUse.Contains(FunctionalGroup))
//                {
//                    if (_CellEnvironment["Realm"][0] == 1.0)
//                    {
//                        // An instance of the terrestrial carbon model class
//                        RevisedTerrestrialPlantModel PlantModel = new RevisedTerrestrialPlantModel();
//
//
//                        // Calculate predicted leaf mass at equilibrium for this stock
//                        double LeafMass = PlantModel.CalculateEquilibriumLeafMass(_CellEnvironment, functionalGroups.GetTraitNames("leaf strategy", FunctionalGroup) == "deciduous");
//
//                        // Initialise the new stock with the relevant properties
//                        NewStock = new Stock((byte)FunctionalGroup, IndividualMass[FunctionalGroup], LeafMass);
//
//                        // Add the new stock to the list of grid cell stocks
//                        _GridCellStocks[FunctionalGroup].Add(NewStock);
//
//                        // Increment the variable tracking the total number of stocks in the model
//                        globalDiagnostics["NumberOfStocksInModel"]++;
//
//
//                    }
//                    else if (FunctionalGroupsToUse.Contains(FunctionalGroup))
//                    {
//                        // Initialise the new stock with the relevant properties
//                        NewStock = new Stock((byte)FunctionalGroup, IndividualMass[FunctionalGroup], 1e12);
//
//                        // Add the new stock to the list of grid cell stocks
//                        _GridCellStocks[FunctionalGroup].Add(NewStock);
//
//                        // Increment the variable tracking the total number of stocks in the model
//                        globalDiagnostics["NumberOfStocksInModel"]++;
//
//                    }
//                    else
//                    {
//                    }
//
//                }
//
//            }
//
//
//
//
//        }
    };
//
//
//}
//
#endif
