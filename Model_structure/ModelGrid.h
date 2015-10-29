#ifndef MODELGRID_H
#define MODELGRID_H
using namespace std;
#include <GridCellCohortHandler.h>
#include <GridCell.h>
#include <map>
#include <limits>
#include <string>

/** \brief A class containing the model grid (composed of individual grid cells) along with grid attributes.
           The model grid is referenced by [Lat index, Lon index]
 */
class ModelGrid {
public:
    //----------------------------------------------------------------------------------------------
    //Variables
    //----------------------------------------------------------------------------------------------
    
    // Variable to make sure that not more than one grid is instantiated
    unsigned NumGrids = 0;

    // Model grid standardised missing value (applied to all grid cells)
    double GlobalMissingValue = -9999;

    /** \brief   minimum latitude of the grid */
    float MinLatitude;
    /** \brief leftmost longitude of the leftmost cell of the grid */
    float MinLongitude;
    /** \brief lowest latitude of the highest cell in the grid */
    float MaxLatitude;
    /** \brief leftmost longitude of the rightmost cell in the grid */
    float MaxLongitude;
    /** \brief latitudinal length of each grid cell. Currently assumes all cells are equal sized. */
    float LatCellSize;
    /** \brief the longitudinal length of each grid cell. Currently assumes all cells are equal sized.  */
    float LonCellSize;
    /** \brief The rarefaction of grid cells to be applied to active cells in the model grid */
    int GridCellRarefaction;
    /** \brief The number of latitudinal cells in the model grid */
    long NumLatCells;
    /** \brief The number of longitudinal cells in the model grid */
    long NumLonCells;
    /** \brief The bottom (southern-most) latitude of each row of grid cells */
    vector<float> Lats;
    /** \brief The left (western-most) longitude of each column of grid cells */
    vector<float> Lons;
    /** \brief Array of grid cells */
    vector< vector <GridCell> > InternalGrid;
    /** \brief
    The heights of grid cells in each latitudinal band
     */
    vector<double> CellHeightsKm;
    /** \brief
    The widths of grid cells in each latitudinal band
     */
    vector<double> CellWidthsKm;
    //
    /** \brief An instance of the simple random number generator class */
    std::default_random_engine RandomNumberGenerator;

    /** \brief Instance of the class to perform general functions */
    UtilityFunctions Utilities;

    //----------------------------------------------------------------------------------------------
    //Methods
    //----------------------------------------------------------------------------------------------
     //empty placeholder constructor just to get compilation to work.

    ModelGrid(){    ;
    }
    //----------------------------------------------------------------------------------------------
    /** \brief Constructor for model grid to construct the grid for specific locations
    @param minLat Minimum grid latitude (degrees) 
    @param minLon Minimum grid longitude (degrees, currently -180 to 180) 
    @param maxLat Maximum grid latitude (degrees) 
    @param maxLon Maximum grid longitude (degrees, currently -180 to 180) 
    @param latCellSize Latitudinal size of grid cells 
    @param lonCellSize Longitudinal size of grid cells 
    @param cellList List of indices of active cells in the model grid 
    @param enviroStack List of environmental data layers 
    @param cohortFunctionalGroups The functional group definitions for cohorts in the model 
    @param stockFunctionalGroups The functional group definitions for stocks in the model 
    @param globalDiagnostics Global diagnostic variables 
    @param nextCohortID The unique ID to assign to the next cohort created 
    @param tracking Whether process tracking is enabled 
     */
    ModelGrid(float minLat, float minLon, float maxLat, float maxLon, float latCellSize, float lonCellSize, vector<vector<unsigned>>&cellList,
            map<string, EnviroData>& enviroStack, FunctionalGroupDefinitions& cohortFunctionalGroups,
            FunctionalGroupDefinitions& stockFunctionalGroups, map<string, double>& globalDiagnostics) {
        // Add one to the counter of the number of grids. If there is more than one model grid, exit the program with a debug crash.
        NumGrids = NumGrids + 1;
        assert(NumGrids < 2 && "You have initialised more than one grid on which to apply models. At present, this is not supported");

        // CURRENTLY DEFINING MODEL CELLS BY BOTTOM LEFT CORNER
        MinLatitude = minLat;
        MinLongitude = minLon;
        MaxLatitude = maxLat;
        MaxLongitude = maxLon;
        LatCellSize = latCellSize;
        LonCellSize = lonCellSize;
        GridCellRarefaction = 1;


        // Check to see if the number of grid cells is an integer
        double u = (MaxLatitude - MinLatitude) / LatCellSize - floor((MaxLatitude - MinLatitude) / LatCellSize);
        assert((u <= 1.e-15) && "Error: number of grid cells is non-integer: check cell size");

        NumLatCells = (long) ((MaxLatitude - MinLatitude) / LatCellSize);
        NumLonCells = (long) ((MaxLongitude - MinLongitude) / LonCellSize);
        Lats.resize(NumLatCells);
        Lons.resize(NumLonCells);

        // Set up latitude and longitude vectors - lower left
        for (int ii = 0; ii < NumLatCells; ii++) {
            Lats[ii] = MinLatitude + ii * LatCellSize;
        }
        for (int jj = 0; jj < NumLonCells; jj++) {
            Lons[jj] = MinLongitude + jj * LonCellSize;
        }

        // Instantiate a grid of grid cells
        InternalGrid.resize(NumLatCells);
        for (auto &g : InternalGrid)g.resize(NumLonCells);

        cout << "Initialising grid cell environment:" << endl;

        int Count = 0;

        int NCells = cellList.size();
   
            // Loop over cells to set up the model grid
            for (int ii = 0; ii < cellList.size(); ii++) {

                // Create the grid cell at the specified position
                InternalGrid[cellList[ii][0]][ cellList[ii][1]] = GridCell(Lats[cellList[ii][0]], cellList[ii][0],
                        Lons[cellList[ii][1]], cellList[ii][1], latCellSize, lonCellSize, enviroStack, GlobalMissingValue,
                        cohortFunctionalGroups, stockFunctionalGroups, globalDiagnostics);

                     Count++;
                    //cout<<"\rInitialised "<<Count<<" of"<< NCells<<endl;

            }
        

        
            //InterpolateMissingValues();//MB data should be prepared so this is not needed

            // Fill in the array of dispersable perimeter lengths for each grid cell
            //CalculatePerimeterLengthsAndCellsDispersableTo();

            CellHeightsKm.resize(Lats.size());
            CellWidthsKm.resize(Lats.size());

            // Calculate the lengths of widths of grid cells in each latitudinal strip
            // Assume that we are at the midpoint of each cell when calculating lengths
            for (int ii = 0; ii < Lats.size(); ii++) {
                CellHeightsKm[ii] = Utilities.CalculateLengthOfDegreeLatitude(Lats[ii] + LatCellSize / 2) * LatCellSize;
                CellWidthsKm[ii] = Utilities.CalculateLengthOfDegreeLongitude(Lats[ii] + LatCellSize / 2) * LonCellSize;
            }
        

        cout << "\n" << endl;
        cout.flush();

    }
   //----------------------------------------------------------------------------------------------
   /** \brief Estimates missing environmental data for grid cells by interpolation */
    void InterpolateMissingValues() {
        map<string, vector<double>> WorkingCellEnvironment;
        bool Changed = false;

        for (unsigned ii = 0; ii < NumLatCells; ii++) {
            for (unsigned jj = 0; jj < NumLonCells; jj++) {
                WorkingCellEnvironment = GetCellEnvironment(ii, jj);

                // If the cell environment does not contain valid NPP data then interpolate values
                if (!InternalGrid[ii][jj].ContainsData(WorkingCellEnvironment["NPP"], WorkingCellEnvironment["Missing Value"][0])) {
                    //If NPP doesn't exist the interpolate from surrounding values (of the same realm)
                    WorkingCellEnvironment["NPP"] = GetInterpolatedValues(ii, jj, GetCellLatitude(ii), GetCellLongitude(jj), "NPP", WorkingCellEnvironment["Realm"][0]);

                    //Calculate NPP seasonality - for use in converting annual NPP estimates to monthly
                    WorkingCellEnvironment["Seasonality"] = InternalGrid[ii][jj].CalculateNPPSeasonality(WorkingCellEnvironment["NPP"], WorkingCellEnvironment["Missing Value"][0]);
                    Changed = true;
                }                    // Otherwise convert the missing data values to zeroes where they exist amongst valid data eg in polar regions.
                else {
                    WorkingCellEnvironment["NPP"] = InternalGrid[ii][jj].ConvertMissingValuesToZero(WorkingCellEnvironment["NPP"], WorkingCellEnvironment["Missing Value"][0]);
                }

                // If the cell environment does not contain valid monthly mean diurnal temperature range data then interpolate values
                if (InternalGrid[ii][jj].ContainsMissingValue(WorkingCellEnvironment["DiurnalTemperatureRange"], WorkingCellEnvironment["Missing Value"][0])) {
                    //If NPP doesn't exist the interpolate from surrounding values (of the same realm)
                    WorkingCellEnvironment["DiurnalTemperatureRange"] = FillWithInterpolatedValues(ii, jj, GetCellLatitude(ii), GetCellLongitude(jj), "DiurnalTemperatureRange", WorkingCellEnvironment["Realm"][0]);

                    Changed = true;
                }

                // Same for u and v velocities
                if (!InternalGrid[ii][jj].ContainsData(WorkingCellEnvironment["uVel"], WorkingCellEnvironment["Missing Value"][0])) {
                    //If u doesn't exist the interpolate from surrounding values (of the same realm)
                    WorkingCellEnvironment["uVel"] = GetInterpolatedValues(ii, jj, GetCellLatitude(ii), GetCellLongitude(jj), "uVel", WorkingCellEnvironment["Realm"][0]);

                    Changed = true;
                }                    // Otherwise convert the missing data values to zeroes where they exist amongst valid data eg in polar regions.
                else {
                    WorkingCellEnvironment["uVel"] = InternalGrid[ii][jj].ConvertMissingValuesToZero(WorkingCellEnvironment["uVel"], WorkingCellEnvironment["Missing Value"][0]);
                }

                if (!InternalGrid[ii][jj].ContainsData(WorkingCellEnvironment["vVel"], WorkingCellEnvironment["Missing Value"][0])) {
                    //If v vel doesn't exist the interpolate from surrounding values (of the same realm)
                    WorkingCellEnvironment["vVel"] = GetInterpolatedValues(ii, jj, GetCellLatitude(ii), GetCellLongitude(jj), "vVel", WorkingCellEnvironment["Realm"][0]);

                    Changed = true;
                }                    // Otherwise convert the missing data values to zeroes where they exist amongst valid data eg in polar regions.
                else {
                    WorkingCellEnvironment["vVel"] = InternalGrid[ii][jj].ConvertMissingValuesToZero(WorkingCellEnvironment["vVel"], WorkingCellEnvironment["Missing Value"][0]);
                }

                if (Changed) InternalGrid[ii][jj].CellEnvironment = WorkingCellEnvironment;
            }
        }
    }
    //----------------------------------------------------------------------------------------------
    /** \brief Calculate the weighted average of surrounding grid cell data, where those grid cells are of the specified realm and contain
    non missing data values
    @param latIndex Index of the latitude cell for which the weighted average over surrounding cells is requested 
    @param lonIndex Index of the longitude cell for which the weighted average over surrounding cells is requested 
    @param lat Latitude of the cell for which the weighted value is requested 
    @param lon Longitude of the cell for which the weighted value is requested 
    @param dataName Names of the data for which weighted value is requested 
    @param realm Realm of the grid cell for which data is to be averaged over 
    @return The weighted average value of the specified data type across surrounding grid cells of the specified realm
     */
    vector<double> GetInterpolatedValues(unsigned latIndex, unsigned lonIndex, double lat, double lon, string dataName, double realm) {
        map<string, vector<double>> TempCellEnvironment = GetCellEnvironment(latIndex, lonIndex);
        vector<double> InterpData(TempCellEnvironment[dataName].size());
        vector<unsigned> InterpCount(TempCellEnvironment[dataName].size());

        unsigned LowerLatIndex = latIndex - 1;
        unsigned UpperLatIndex = latIndex + 1;
        unsigned LowerLonIndex = lonIndex - 1;
        unsigned UpperLonIndex = lonIndex + 1;


        if (latIndex == 0) LowerLatIndex = latIndex;
        if (lat == MaxLatitude) UpperLatIndex = latIndex;

        if (lonIndex == 0) LowerLonIndex = lonIndex;
        if (lon == MaxLongitude) UpperLonIndex = lonIndex;

        //Loop over surrounding cells in the datalayer
        for (unsigned ii = LowerLatIndex; ii <= UpperLatIndex; ii++) {
            for (unsigned jj = LowerLonIndex; jj < UpperLonIndex; jj++) {
                if (ii < NumLatCells && jj < NumLonCells) {
                    TempCellEnvironment = GetCellEnvironment(ii, jj);

                    for (unsigned hh = 0; hh < InterpData.size(); hh++) {
                        //If the cell contains data then sum this and increment count
                        if (TempCellEnvironment[dataName][hh] != TempCellEnvironment["Missing Value"][0] && TempCellEnvironment["Realm"][0] == realm) {
                            InterpData[hh] += TempCellEnvironment[dataName][hh];
                            InterpCount[hh]++;
                        }
                    }
                }
            }
        }

        //take the mean over surrounding valid cells for each timestep
        for (int hh = 0; hh < InterpData.size(); hh++) {
            if (InterpCount[hh] > 0) {
                InterpData[hh] /= InterpCount[hh];
            } else {
                InterpData[hh] = 0.0;
            }
        }
        return InterpData;
    }
    //----------------------------------------------------------------------------------------------
    /** \brief Calculate the weighted average of surrounding grid cell data, where those grid cells are of the specified realm and contain
    non missing data values
    @param latIndex Index of the latitude cell for which the weighted average over surrounding cells is requested 
    @param lonIndex Index of the longitude cell for which the weighted average over surrounding cells is requested 
    @param lat Latitude of the cell for which the weighted value is requested 
    @param lon Longitude of the cell for which the weighted value is requested 
    @param dataName Names of the data for which weighted value is requested 
    @param realm Realm of the grid cell for which data is to be averaged over 
    @returns The weighted average value of the specified data type across surrounding grid cells of the specified realm
     */
    vector<double> FillWithInterpolatedValues(unsigned latIndex, unsigned lonIndex, double lat, double lon, string dataName, double realm) {
        map<string, vector<double> > TempCellEnvironment = GetCellEnvironment(latIndex, lonIndex);
        vector<double> InterpData(TempCellEnvironment[dataName].size());
        vector<unsigned> InterpCount(TempCellEnvironment[dataName].size());
        unsigned LowerLatIndex = latIndex - 1;
        unsigned UpperLatIndex = latIndex + 1;
        unsigned LowerLonIndex = lonIndex - 1;
        unsigned UpperLonIndex = lonIndex + 1;


        if (latIndex == 0) LowerLatIndex = latIndex;
        if (lat == MaxLatitude) UpperLatIndex = latIndex;

        if (lonIndex == 0) LowerLonIndex = lonIndex;
        if (lon == MaxLongitude) UpperLonIndex = lonIndex;

        for (unsigned hh = 0; hh < InterpData.size(); hh++) {
            if (TempCellEnvironment[dataName][hh] == TempCellEnvironment["Missing Value"][0]) {
                //Loop over surrounding cells in the datalayer
                for (unsigned ii = LowerLatIndex; ii <= UpperLatIndex; ii++) {
                    for (unsigned jj = LowerLonIndex; jj <= UpperLonIndex; jj++) {
                        if (ii < NumLatCells && jj < NumLonCells) {
                            TempCellEnvironment = GetCellEnvironment(ii, jj);

                            //If the cell contains data then sum this and increment count
                            if (TempCellEnvironment[dataName][hh] != TempCellEnvironment["Missing Value"][0] && TempCellEnvironment["Realm"][0] == realm) {
                                InterpData[hh] += TempCellEnvironment[dataName][hh];
                                InterpCount[hh]++;
                            }

                        }
                    }
                }
                //take the mean over surrounding valid cells for each timestep
                if (InterpCount[hh] > 0) {
                    InterpData[hh] /= InterpCount[hh];
                } else {
                    InterpData[hh] = 0.0;
                }
            } else {
                InterpData[hh] = TempCellEnvironment[dataName][hh];
            }
        }


        return InterpData;
    }
    //----------------------------------------------------------------------------------------------
    /** \brief Returns the stocks within the specified grid cell
    @param latIndex Latitude index 
    @param lonIndex Longitude index 
    @returns The stock handler for the specified grid cell
     */
    GridCellStockHandler& GetGridCellStocks(unsigned latIndex, unsigned lonIndex) {
        return InternalGrid[latIndex][ lonIndex].GridCellStocks;
    }
    //----------------------------------------------------------------------------------------------
    /** \brief Sets the stocks in the specified grid cell to the passed stocks
    @param newGridCellStocks New stocks for the grid cell 
    @param latIndex Latitude index 
    @param lonIndex Longitude index 
     */
    void SetGridCellStocks(GridCellStockHandler newGridCellStocks, unsigned latIndex, unsigned lonIndex) {
        InternalGrid[latIndex][ lonIndex].GridCellStocks = newGridCellStocks;
    }
    //----------------------------------------------------------------------------------------------
    /** \brief
    Returns the array (indexed by functional group) of lists of gridCellCohorts for the specified grid cell
    @param latIndex Latitude index of grid cell 
    @param lonIndex Longitude index of grid cell 
    @returns Array (indexed by functional group) of lists of gridCellCohorts
     */
    GridCellCohortHandler& GetGridCellCohorts(unsigned latIndex, unsigned lonIndex) {
        return InternalGrid[latIndex][ lonIndex].GridCellCohorts;
    }
    //----------------------------------------------------------------------------------------------
    /** \brief
    Extracts an individual cohort from a particular grid cell
    @param latIndex Latitude index of grid cell 
    @param lonIndex Longitude index of grid cell 
    @param functionalGroup Functional group of cohort 
    @param positionInList Index of cohort position in the list 
    @return what?
     */
    Cohort& GetGridCellIndividualCohort(unsigned latIndex, unsigned lonIndex, int functionalGroup, int positionInList) {
        return InternalGrid[latIndex][ lonIndex].GridCellCohorts[functionalGroup].at(positionInList);
    }
    //----------------------------------------------------------------------------------------------
    /** \brief Remove an individual cohort from a functional group; necessary due to dispersal moving cohorts from one cell to another
    @param functionalGroup Cohort functional group 
     */
    void DeleteGridCellIndividualCohort(Cohort c) {
      
        auto begin = InternalGrid[c.origin[0]][c.origin[1]].GridCellCohorts[c.FunctionalGroupIndex].begin();
        InternalGrid[c.origin[0]][c.origin[1]].GridCellCohorts[c.FunctionalGroupIndex].erase(begin + c.positionInList);

    }
    //----------------------------------------------------------------------------------------------
    /** \brief Replace the gridCellCohorts in a grid cell with a new list of gridCellCohorts
    @param newGridCellCohorts The new list of gridCellCohorts 
    @param latIndex Grid cell latitude index 
    @param lonIndex Grid cell longitude index 
     */
    void SetGridCellCohorts(GridCellCohortHandler newGridCellCohorts, unsigned latIndex, unsigned lonIndex) {
        InternalGrid[latIndex][ lonIndex].GridCellCohorts = newGridCellCohorts;
    }
    //----------------------------------------------------------------------------------------------
    /** \brief Add a new cohort to an existing list of cohorts in the grid cell - or create a new list if there is not one present
    @param latIndex Latitude index of the grid cell 
    @param lonIndex Longitude index of the grid cell 
    @param functionalGroup Functional group of the cohort (i.e. array index) 
    @param cohortToAdd The cohort object to add 
     */
    void AddNewCohortToGridCell(Cohort c) {
        c.origin[0]=c.destination[0];c.origin[1]=c.destination[1];
        c.positionInList=InternalGrid[c.destination[0]][ c.destination[1]].GridCellCohorts[c.FunctionalGroupIndex].size();
        InternalGrid[c.destination[0]][ c.destination[1]].GridCellCohorts[c.FunctionalGroupIndex].push_back(c);
    }
    //----------------------------------------------------------------------------------------------
    /** \brief Return the value of a specified environmental layer from an individual grid cell    @param variableName The name of the environmental lyaer 
    @param timeInterval The desired time interval within the environmental variable (i.e. 0 if it is a yearly variable
    or the month index - 0=Jan, 1=Feb etc. - for monthly variables) 
    @param latCellIndex The latitudinal cell index 
    @param lonCellIndex The longitudinal cell index 
    @param variableExists Returns false if the environmental layer does not exist, true if it does 
    @returns The value of the environmental layer, or a missing value if the environmental layer does not exist
     */
    //MB NB returns a COPY of the value stored

    double GetEnviroLayer(string variableName, unsigned timeInterval, unsigned latCellIndex, unsigned lonCellIndex, bool variableExists) {
        return InternalGrid[latCellIndex][ lonCellIndex].GetEnviroLayer(variableName, timeInterval, variableExists);
    }
    //----------------------------------------------------------------------------------------------
    /** \brief
    Set the value of a specified environmental layer in an individual grid cell

    @param variableName The name of the environmental layer 
    @param timeInterval The time interval within the environmental variable to set (i.e. 0 if it is a yearly variable
    or the month index - 0=Jan, 1=Feb etc. - for monthly variables) 
    @param setValue The value to set 
    @param latCellIndex The latitudinal cell index 
    @param lonCellIndex The longitudinal cell index 
    @return True if the value is set successfully, false otherwise
     */
    bool SetEnviroLayer(string variableName, unsigned timeInterval, double setValue, unsigned latCellIndex, unsigned lonCellIndex) {
        return InternalGrid[latCellIndex][ lonCellIndex].SetEnviroLayer(variableName, timeInterval, setValue);
    }
    //----------------------------------------------------------------------------------------------
    /** \brief Set the value of a given delta type for the specified ecological process within the specified grid cell
    @param deltaType The type of delta value to set (e.g. 'biomass', 'abundance' etc.) 
    @param ecologicalProcess The name of the ecological process to set the value of delta for 
    @param setValue The value to set 
    @param latCellIndex The latitudinal index of the cell 
    @param lonCellIndex The longitudinal index of the cell 
    @returns True if the value is set successfully, false otherwise
     */
    bool SetDeltas(string deltaType, string ecologicalProcess, double setValue, unsigned latCellIndex, unsigned lonCellIndex) {
        return InternalGrid[latCellIndex][ lonCellIndex].SetDelta(deltaType, ecologicalProcess, setValue);
    }
    //----------------------------------------------------------------------------------------------
    /** \brief
    Get the total of a state variable for specific cells
    @param variableName The name of the variable 
    @param functionalGroups A vector of functional group indices to consider 
    @param cellIndices List of indices of active cells in the model grid 
    @param stateVariableType A string indicating the type of state variable; 'cohort' or 'stock' 
    @returns Summed value of variable over whole grid
     */
    double StateVariableGridTotal(string variableName, string traitValue, vector<int> functionalGroups, vector<vector<unsigned>> cellIndices, string stateVariableType, MadingleyModelInitialisation initialisation) {

        double tempVal = 0;

        vector<vector< double> > TempStateVariable = GetStateVariableGrid(variableName, traitValue, functionalGroups, cellIndices, stateVariableType, initialisation);

        // Loop through and sum values across a grid, excluding missing values
        for (int ii = 0; ii < cellIndices.size(); ii++) {
            tempVal += TempStateVariable[cellIndices[ii][0]][ cellIndices[ii][1]];
        }

        return tempVal;
    }
    //----------------------------------------------------------------------------------------------
    /** \brief Gets a state variable for specified functional groups of specified entity types in a specified grid cell
    @param variableName The name of the variable to get: 'biomass' or 'abundance' 
    @param functionalGroups The functional group indices to get the state variable for 
    @param latCellIndex The latitudinal index of the cell 
    @param lonCellIndex The longitudinal index of the cell 
    @param stateVariableType The type of entity to return the state variable for: 'stock' or 'cohort' 
    @return The state variable for specified functional groups of specified entity types in a specified grid cell
     */
    double GetStateVariable(string variableName, string traitValue, vector<int> functionalGroups, unsigned latCellIndex, unsigned lonCellIndex, string stateVariableType, MadingleyModelInitialisation modelInitialisation) {

        double returnValue = 0.0;
        map<string, int> vn, sv;
        sv["cohort"] = 0;
        sv["stock"] = 1;
        vn["biomass" ] = 0;
        vn["abundance"] = 1;

        //lowercase the string - a bit clunky...but then C++ strings are a bit
        transform(variableName.begin(), variableName.end(), variableName.begin(), ::tolower);
        transform(stateVariableType.begin(), stateVariableType.end(), stateVariableType.begin(), ::tolower);
        GridCellCohortHandler TempCohorts = InternalGrid[latCellIndex][ lonCellIndex].GridCellCohorts;
        GridCellStockHandler TempStocks = InternalGrid[latCellIndex][ lonCellIndex].GridCellStocks;

        switch (sv[stateVariableType]) {
            case 0://"cohort":


                switch (vn[variableName]) {
                    case 0://"biomass":
                        if (traitValue != "Zooplankton") {
                            for (int f : functionalGroups) {
                                for (auto item : TempCohorts[f]) {
                                    returnValue += ((item.IndividualBodyMass + item.IndividualReproductivePotentialMass) * item.CohortAbundance);
                                }
                            }
                        } else {
                            for (int f : functionalGroups) {
                                for (auto item : TempCohorts[f]) {
                                    if (item.IndividualBodyMass <= modelInitialisation.PlanktonDispersalThreshold)
                                        returnValue += ((item.IndividualBodyMass + item.IndividualReproductivePotentialMass) * item.CohortAbundance);
                                }
                            }
                        }
                        break;

                    case 1://"abundance":
                        if (traitValue != "Zooplankton") {
                            for (int f : functionalGroups) {
                                for (auto item : TempCohorts[f]) {
                                    returnValue += item.CohortAbundance;
                                }
                            }
                        } else {
                            for (int f : functionalGroups) {
                                for (auto item : TempCohorts[f]) {
                                    if (item.IndividualBodyMass <= modelInitialisation.PlanktonDispersalThreshold)
                                        returnValue += item.CohortAbundance;
                                }
                            }
                        }
                        break;

                    default:
                        cout << "For cohorts, state variable name must be either 'biomass' or 'abundance'";
                        exit(1);
                        break;
                }
                break;

            case 1://"stock":

                switch (vn[variableName]) {
                    case 0://"biomass":
                        for (int f : functionalGroups) {
                            for (auto item : TempStocks[f]) {
                                returnValue += item.TotalBiomass;
                            }
                        }
                        break;
                    default:
                        cout << "For stocks, state variable name must be 'biomass'" << endl;
                        exit(1);
                        break;
                }
                break;

            default:
                cout << "State variable type must be either 'cohort' or 'stock'" << endl;
                exit(1);
                break;

        }
        return returnValue;
    }
    //----------------------------------------------------------------------------------------------
    /** \brief Gets a state variable density for specified functional groups of specified entity types in a specified grid cell
    @param variableName The name of the variable to get: 'biomass' or 'abundance' 
    @param functionalGroups The functional group indices to get the state variable for 
    @param latCellIndex The latitudinal index of the cell 
    @param lonCellIndex The longitudinal index of the cell 
    @param stateVariableType The type of entity to return the state variable for: 'stock' or 'cohort' 
    @return The state variable density for specified functional groups of specified entity types in a specified grid cell
     */
    double GetStateVariableDensity(string variableName, string traitValue, vector<int> functionalGroups, unsigned latCellIndex, unsigned lonCellIndex, string stateVariableType, MadingleyModelInitialisation modelInitialisation) {

        double returnValue = 0.0;
        map<string, int> vn, sv;
        sv["cohort"] = 0;
        sv["stock"] = 1;
        vn["biomass" ] = 0;
        vn["abundance"] = 1;

        //lowercase the string - a bit clunky...but then C++ strings are a bit
        transform(variableName.begin(), variableName.end(), variableName.begin(), ::tolower);
        transform(stateVariableType.begin(), stateVariableType.end(), stateVariableType.begin(), ::tolower);

        GridCellStockHandler TempStocks = InternalGrid[latCellIndex][ lonCellIndex].GridCellStocks;
        GridCellCohortHandler TempCohorts = InternalGrid[latCellIndex][ lonCellIndex].GridCellCohorts;

        switch (sv[stateVariableType]) {
            case 0://"cohort":


                switch (vn[variableName]) {
                    case 0://"biomass":
                        if (traitValue != "Zooplankton (all)") {
                            for (int f : functionalGroups) {
                                for (auto item : TempCohorts[f]) {
                                    returnValue += ((item.IndividualBodyMass + item.IndividualReproductivePotentialMass) * item.CohortAbundance);
                                }
                            }
                        } else {
                            for (int f : functionalGroups) {
                                for (auto item : TempCohorts[f]) {
                                    if (item.IndividualBodyMass <= modelInitialisation.PlanktonDispersalThreshold)
                                        returnValue += ((item.IndividualBodyMass + item.IndividualReproductivePotentialMass) * item.CohortAbundance);
                                }
                            }
                        }
                        break;

                    case 1://"abundance":
                        if (traitValue != "Zooplankton (all)") {
                            for (int f : functionalGroups) {
                                for (auto item : TempCohorts[f]) {
                                    returnValue += item.CohortAbundance;
                                }
                            }
                        } else {
                            for (int f : functionalGroups) {
                                for (auto item : TempCohorts[f]) {
                                    if (item.IndividualBodyMass <= modelInitialisation.PlanktonDispersalThreshold)
                                        returnValue += item.CohortAbundance;
                                }
                            }
                        }
                        break;

                    default:
                        cout << "For cohorts, state variable name must be either 'biomass' or 'abundance'" << endl;
                        exit(1);
                        break;
                }
                break;

            case 1://"stock":

                switch (vn[variableName]) {
                    case 0://"biomass":
                        for (int f : functionalGroups) {
                            for (auto item : TempStocks[f]) {
                                returnValue += item.TotalBiomass;
                            }
                        }
                        break;
                    default:
                        cout << "For stocks, state variable name must be 'biomass'" << endl;
                        exit(1);
                        break;
                }
                break;

            default:
                cout << "State variable type must be either 'cohort' or 'stock'" << endl;
                exit(1);
                break;

        }


        return returnValue / (InternalGrid[latCellIndex][ lonCellIndex].CellEnvironment["Cell Area"][0]);
    }
    //----------------------------------------------------------------------------------------------
    /** \brief Get the mean density of a state variable for specific cells
    @param variableName The name of the variable 
    @param functionalGroups A vector of functional group indices to consider 
    @param cellIndices List of indices of active cells in the model grid 
    @param stateVariableType A string indicating the type of state variable; 'cohort' or 'stock' 
    @return Mean density of variable over whole grid
     */
    double StateVariableGridMeanDensity(string variableName, string traitValue, vector<int> functionalGroups, vector<vector<unsigned>> cellIndices, string stateVariableType, MadingleyModelInitialisation initialisation) {

        double tempVal = 0;

        vector<vector< double> > TempStateVariable = GetStateVariableGridDensityPerSqKm(variableName, traitValue, functionalGroups, cellIndices, stateVariableType, initialisation);

        // Loop through and sum values across a grid, excluding missing values
        for (int ii = 0; ii < cellIndices.size(); ii++) {
            tempVal += TempStateVariable[cellIndices[ii][0]][ cellIndices[ii][1]];
        }

        return tempVal / cellIndices.size();
    }
    //----------------------------------------------------------------------------------------------
    /** \brief Return an array of values for a single state variable over specific cells
    @param variableName Variable name 
    @param functionalGroups A vector of functional group indices to consider 
    @param cellIndices List of indices of active cells in the model grid 
    @param stateVariableType A string indicating the type of state variable; 'cohort' or 'stock' 
    @return Array of state variable values for each grid cell
     */
    //MB NB ***!!! array bounds checking not yet in place
    vector<vector< double> > GetStateVariableGrid(string variableName, string traitValue, vector<int> functionalGroups, vector<vector<unsigned>> cellIndices, string stateVariableType, MadingleyModelInitialisation initialisation) {
        vector<vector< double> > TempStateVariable(NumLatCells);
        for (auto &t : TempStateVariable)t.resize(NumLonCells);
        map<string, int> vn;
        vn["biomass" ] = 0;
        vn["abundance"] = 1;

        //lowercase the string - a bit clunky...but then C++ strings are a bit
        transform(variableName.begin(), variableName.end(), variableName.begin(), ::tolower);
        transform(stateVariableType.begin(), stateVariableType.end(), stateVariableType.begin(), ::tolower);
        switch (vn[variableName]) {
            case 0://"biomass":
                for (int ii = 0; ii < cellIndices.size(); ii++) {
                    // Check whether the state variable concerns cohorts or stocks
                    if (stateVariableType == "cohort") {
                        if (traitValue != "Zooplankton") {
                            //                                // Check to make sure that the cell has at least one cohort
                            //                                if (InternalGrid[cellIndices[ii][0]][ cellIndices[ii][1]].GridCellCohorts != null)
                            //                                {
                            for (int nn = 0; nn < functionalGroups.size(); nn++) {
                                //                                        if (InternalGrid[cellIndices[ii][0]][ cellIndices[ii][1]].GridCellCohorts[functionalGroups[nn]] != null)
                                //                                        {
                                for (Cohort item : InternalGrid[cellIndices[ii][0]][ cellIndices[ii][1]].GridCellCohorts[functionalGroups[nn]]) {
                                    TempStateVariable[cellIndices[ii][0]][ cellIndices[ii][1]] += ((item.IndividualBodyMass + item.IndividualReproductivePotentialMass) * item.CohortAbundance);
                                }
                                //                                        }
                            }
                            //                                }
                        } else {
                            //                                //Check to make sure that the cell has at least one cohort
                            //                                if (InternalGrid[cellIndices[ii][0], cellIndices[ii][1]].GridCellCohorts != null)
                            //                                {
                            for (int nn = 0; nn < functionalGroups.size(); nn++) {
                                //                                        if (InternalGrid[cellIndices[ii][0]][ cellIndices[ii][1]].GridCellCohorts[functionalGroups[nn]] != null)
                                //                                        {
                                for (Cohort item : InternalGrid[cellIndices[ii][0]][ cellIndices[ii][1]].GridCellCohorts[functionalGroups[nn]]) {
                                    if (item.IndividualBodyMass <= initialisation.PlanktonDispersalThreshold)
                                        TempStateVariable[cellIndices[ii][0]][ cellIndices[ii][1]] += ((item.IndividualBodyMass + item.IndividualReproductivePotentialMass) * item.CohortAbundance);
                                }
                                //                                        }
                            }
                            //                                }
                        }
                    } else if (stateVariableType == "stock") {
                        //                            // Check to make sure that the cell has at least one stock
                        //                            if (InternalGrid[cellIndices[ii][0]][ cellIndices[ii][1]].GridCellStocks != null)
                        //                            {
                        for (int nn = 0; nn < functionalGroups.size(); nn++) {
                            //                                    if (InternalGrid[cellIndices[ii][0]][ cellIndices[ii][1]].GridCellStocks[functionalGroups[nn]] != null)
                            //                                    {
                            for (Stock item : InternalGrid[cellIndices[ii][0]][ cellIndices[ii][1]].GridCellStocks[functionalGroups[nn]]) {
                                TempStateVariable[cellIndices[ii][0]][ cellIndices[ii][1]] += (item.TotalBiomass);

                            }
                            //                                    }
                            // 
                        }
                        //                            }
                    } else {
                        cout << "Variable 'state variable type' must be either 'stock' 'or 'cohort'" << endl;
                        exit(1);
                    }

                }
                break;
            case 1://"abundance":
                for (int ii = 0; ii < cellIndices.size(); ii++) {
                    // Check whether the state variable concerns cohorts or stocks
                    if (stateVariableType == "cohort") {
                        if (traitValue != "Zooplankton") {
                            //                                // Check to make sure that the cell has at least one cohort
                            //                                if (InternalGrid[cellIndices[ii][0], cellIndices[ii][1]].GridCellCohorts != null)
                            //                                {
                            for (int nn = 0; nn < functionalGroups.size(); nn++) {
                                //                                        if (InternalGrid[cellIndices[ii][0]][ cellIndices[ii][1]].GridCellCohorts[functionalGroups[nn]] != null)
                                //                                        {
                                for (Cohort item : InternalGrid[cellIndices[ii][0]][ cellIndices[ii][1]].GridCellCohorts[functionalGroups[nn]]) {
                                    TempStateVariable[cellIndices[ii][0]][ cellIndices[ii][1]] += item.CohortAbundance;
                                }
                                //                                        }
                            }
                            //                                }
                        } else {
                            //                                // Check to make sure that the cell has at least one cohort
                            //                                if (InternalGrid[cellIndices[ii][0], cellIndices[ii][1]].GridCellCohorts != null)
                            //                                {
                            for (int nn = 0; nn < functionalGroups.size(); nn++) {
                                //                                        if (InternalGrid[cellIndices[ii][0]][ cellIndices[ii][1]].GridCellCohorts[functionalGroups[nn]] != null)
                                //                                        {
                                for (Cohort item : InternalGrid[cellIndices[ii][0]][ cellIndices[ii][1]].GridCellCohorts[functionalGroups[nn]]) {
                                    if (item.IndividualBodyMass <= initialisation.PlanktonDispersalThreshold)
                                        TempStateVariable[cellIndices[ii][0]][ cellIndices[ii][1]] += item.CohortAbundance;
                                }
                                //                                        }
                            }
                            //                                }
                        }
                    } else {
                        cout << "Currently abundance cannot be calculated for grid cell stocks" << endl;
                        exit(1);
                    }
                }
                break;
            default:
                cout << "Invalid search string passed for cohort property" << endl;
                exit(1);
                break;
        }
        //
        //            return TempStateVariable;

    }
    //----------------------------------------------------------------------------------------------
    /** \brief  Return an array of values for a single state variable over specific cells, given in densities per km^2
    @param variableName Variable name 
    @param functionalGroups A vector of functional group indices to consider 
    @param cellIndices List of indices of active cells in the model grid 
    @param stateVariableType A string indicating the type of state variable; 'cohort' or 'stock' 
    @return Array of state variable values for each grid cell
     */
    vector<vector< double> > GetStateVariableGridDensityPerSqKm(string variableName, string traitValue, vector<int> functionalGroups, vector<vector<unsigned>> cellIndices, string stateVariableType, MadingleyModelInitialisation initialisation) {
        vector<vector< double> > TempStateVariable(NumLatCells);
        for (auto &t : TempStateVariable)t.resize(NumLonCells);
        double CellArea;

        TempStateVariable = GetStateVariableGrid(variableName, traitValue, functionalGroups, cellIndices, stateVariableType, initialisation);

        for (int ii = 0; ii < cellIndices.size(); ii++) {
            CellArea = GetCellEnvironment(cellIndices[ii][0], cellIndices[ii][1])["Cell Area"][0];
            TempStateVariable[cellIndices[ii][0]][ cellIndices[ii][1]] /= CellArea;
        }

        return TempStateVariable;
    }
    //----------------------------------------------------------------------------------------------
    /** \brief   Return an array of log(values + 1) for a state variable for particular functional groups over specific cells. State variable (currently only biomass or abundance) must be >= 0 in all grid cells
    @param variableName The name of the variable 
    @param functionalGroups A vector of functional group indices to consider 
    @param cellIndices List of indices of active cells in the model grid 
    @param stateVariableType A string indicating the type of state variable; 'cohort' or 'stock' 
    @return Array of log(state variable values +1 ) for each grid cell
     */
    vector<vector< double> > GetStateVariableGridLog(string variableName, string traitValue, vector<int> functionalGroups, vector<vector<unsigned>> cellIndices, string stateVariableType, MadingleyModelInitialisation initialisation) {

        vector<vector< double> > TempStateVariable(NumLatCells);
        for (auto &t : TempStateVariable)t.resize(NumLonCells);

        TempStateVariable = GetStateVariableGrid(variableName, traitValue, functionalGroups, cellIndices, stateVariableType, initialisation);

        for (int ii = 0; ii < cellIndices.size(); ii++) {
            TempStateVariable[cellIndices[ii][0]][ cellIndices[ii][1]] = log(TempStateVariable[cellIndices[ii][0]][ cellIndices[ii][1]] + 1);
        }

        return TempStateVariable;
    }
    //----------------------------------------------------------------------------------------------
    /** \brief  Return an array of log(values + 1) for a state variable for particular functional groups over specific cells. State variable (currently only biomass or abundance) must be >= 0 in all grid cells
    @param variableName The name of the variable 
    @param functionalGroups A vector of functional group indices to consider 
    @param cellIndices List of indices of active cells in the model grid 
    @param stateVariableType A string indicating the type of state variable; 'cohort' or 'stock' 
    @return Array of log(state variable values +1 ) for each grid cell
     */
    vector<vector< double> > GetStateVariableGridLogDensityPerSqKm(string variableName, string traitValue, vector<int> functionalGroups, vector<vector<unsigned>> cellIndices, string stateVariableType, MadingleyModelInitialisation initialisation) {

        vector<vector< double> > TempStateVariable(NumLatCells);
        for (auto &t : TempStateVariable)t.resize(NumLonCells);
        double CellArea;

        TempStateVariable = GetStateVariableGrid(variableName, traitValue, functionalGroups, cellIndices, stateVariableType, initialisation);

        for (int ii = 0; ii < cellIndices.size(); ii++) {
            CellArea = GetCellEnvironment(cellIndices[ii][0], cellIndices[ii][1])["Cell Area"][0];
            TempStateVariable[cellIndices[ii][0]][ cellIndices[ii][1]] /= CellArea;
            TempStateVariable[cellIndices[ii][0]][ cellIndices[ii][1]] = log(TempStateVariable[cellIndices[ii][0]][ cellIndices[ii][1]] + 1);
        }

        return TempStateVariable;

    }
    //----------------------------------------------------------------------------------------------
    /** \brief  Returns, for a given longitude, the appropriate longitude index in the grid
    ASSUMES THAT LONGITUDES IN THE MODEL GRID OBJECT REFER TO LOWER LEFT CORNERS!!!
    @param myLon Longitude, in degrees 
    @return longitude index in the model grid
     */
    unsigned GetLonIndex(double myLon) {
        assert((myLon >= MinLongitude && myLon < MaxLongitude) && "Error: latitude out of range");

        return (unsigned) floor((myLon - MinLongitude) / LonCellSize);

    }
    //----------------------------------------------------------------------------------------------
    /** \brief  Return the longitude of a cell at a particular lon. index
    @param cellLonIndex The longitudinal index (i.e. row) of the cell 
    @return Returns the longitude of the bottom of the cell, in degrees
     */
    double GetCellLongitude(unsigned cellLonIndex) {
        assert((cellLonIndex <= (NumLonCells - 1)) && "Error: Cell index out of range when trying to find the longitude for a particular cell");

        double TempLongitude = std::numeric_limits<unsigned>::max();
        ;

        for (int ii = 0; ii < NumLatCells; ii++) {
            if (InternalGrid.size() >= ii && InternalGrid[cellLonIndex].size() >= cellLonIndex)
                TempLongitude = InternalGrid[ii][ cellLonIndex].Longitude;
        }

        assert(TempLongitude != std::numeric_limits<unsigned>::max() && "Error trying to find cell longitude - no grid cells have been initialised for this latitude index: ");

        return TempLongitude;
    }
    //----------------------------------------------------------------------------------------------
    /** \brief Return the latitude of a cell at a particular lat. index
    @param cellLatIndex The latitudinal index (i.e. row) of the cell 
    @return Returns the latitude of the bottom of the cell, in degrees
     */
    double GetCellLatitude(unsigned cellLatIndex) {
        assert((cellLatIndex <= (NumLatCells - 1)) && "Error: Cell index out of range when trying to find the latitude for a particular cell");

        double TempLatitude = std::numeric_limits<unsigned>::max();

        for (int jj = 0; jj < NumLonCells; jj++) {
            if (InternalGrid.size() >= cellLatIndex && InternalGrid[cellLatIndex].size() >= jj) {
                TempLatitude = InternalGrid[cellLatIndex][ jj].Latitude;
                break;
            }
        }

        assert((TempLatitude != std::numeric_limits<unsigned>::max()) && "Error trying to find cell latitude - no grid cells have been initialised for this latitude index: ");

        return TempLatitude;

    }
    //----------------------------------------------------------------------------------------------
    /** \brief   Returns, for a given latitude, the appropriate latitude index in the grid
    ASSUMES THAT LATITUDES IN THE MODEL GRID OBJECT REFER TO LOWER LEFT CORNERS!!!
    @param myLat Latitude, in degrees 
    @return latitude index in the model grid
     */
    unsigned GetLatIndex(double myLat) {

        assert((myLat >= MinLatitude && myLat < MaxLatitude) && "Error: latitude out of range");

        return (unsigned) floor((myLat - MinLatitude) / LatCellSize);

    }
    //----------------------------------------------------------------------------------------------
    /** \brief   A method to return the values for all environmental data layers for a particular grid cell
    @param cellLatIndex Latitude index of grid cell 
    @param cellLonIndex Longitude index of grid cell 
    @return A sorted list containing environmental data layer names and values
     */
    map<string, vector<double> >& GetCellEnvironment(unsigned cellLatIndex, unsigned cellLonIndex) {
        return InternalGrid[cellLatIndex][ cellLonIndex].CellEnvironment;
    }
    //----------------------------------------------------------------------------------------------
    /** \brief
    A method to return delta values for the specified delta type in a particular grid cell
    @param deltaType The delta type to return 
    @param cellLatIndex Latitude index of grid cell 
    @param cellLonIndex Longitude index of grid cell 
    @return A sorted list containing deltas
     */
    map<string, double>& GetCellDeltas(string deltaType, unsigned cellLatIndex, unsigned cellLonIndex) {
        return InternalGrid[cellLatIndex][ cellLonIndex].Deltas[deltaType];
    }
    //----------------------------------------------------------------------------------------------
    /** \brief A method to return all delta values in a particular grid cell
    @param cellLatIndex Latitude index of grid cell 
    @param cellLonIndex Longitude index of grid cell 
    @return A sorted list of sorted lists containing deltas
     */
    map<string, map<string, double>>&GetCellDeltas(unsigned cellLatIndex, unsigned cellLonIndex) {
        return InternalGrid[cellLatIndex][ cellLonIndex].Deltas;
    }
    //----------------------------------------------------------------------------------------------
    /** \brief  Get a grid of values for an environmental data layer
    @param enviroVariable  The name of the environmental data layer 
    @param timeInterval The desired time interval within the environmental variable (i.e. 0 if it is a yearly variable
    or the month index - 0=Jan, 1=Feb etc. - for monthly variables) 
    @return The values in each grid cell
     */
    vector<vector< double> > GetEnviroGrid(string enviroVariable, unsigned timeInterval) {
        // Check to see if environmental variable exists
        for (int ii = 0; ii < NumLatCells; ii++) {
            for (int jj = 0; jj < NumLonCells; jj++) {
                if (InternalGrid.size() >= ii && InternalGrid[ii].size() >= jj)
                    assert(InternalGrid[ii][ jj].CellEnvironment.count(enviroVariable) != 0 && "Environmental variable not found when running GetEnviroGrid");
            }
        }

        vector<vector< double> > outputData(NumLatCells);
        for (auto e : outputData)e.resize(NumLonCells);

        for (int ii = 0; ii < NumLatCells; ii += GridCellRarefaction) {
            for (int jj = 0; jj < NumLonCells; jj += GridCellRarefaction) {
                outputData[ii][ jj] = InternalGrid[ii][ jj].CellEnvironment[enviroVariable][timeInterval];
            }
        }

        return outputData;
    }
    //----------------------------------------------------------------------------------------------
    /** \brief
    Get a grid of values for an environmental data layer in specific cells
    @param enviroVariable The name of the environmental data layer to return 
    @param timeInterval The desired time interval for which to get data (i.e. 0 if it is a yearly variable
    or the month index - 0=Jan, 1=Feb etc. - for monthly variables) 
    @param cellIndices List of active cells in the model grid 
    @return The values in each grid cell
     */
    vector<vector< double> > GetEnviroGrid(string enviroVariable, unsigned timeInterval, vector<vector<unsigned>> cellIndices) {
        // Check to see if environmental variable exists
        for (int ii = 0; ii < cellIndices.size(); ii++) {
            if (InternalGrid.size() >= cellIndices[ii][0] && InternalGrid[cellIndices[ii][0]].size() >= cellIndices[ii][1])
                assert(InternalGrid[cellIndices[ii][0]][ cellIndices[ii][1]].CellEnvironment.count(enviroVariable) != 0
                    && "Environmental variable not found when running GetEnviroGrid");
        }

        // Create grid to hold the data to return
        vector<vector< double> > outputData(NumLatCells);
        for (auto e : outputData)e.resize(NumLonCells);

        for (int ii = 0; ii < cellIndices.size(); ii++) {
            outputData[cellIndices[ii][0]][ cellIndices[ii][1]] = InternalGrid[cellIndices[ii][0]][ cellIndices[ii][1]].CellEnvironment
                    [enviroVariable][timeInterval];
        }

        return outputData;
    }
    //----------------------------------------------------------------------------------------------
    /** \brief  Return the total over the whole grid for an environmental variable
    @param enviroVariable The environmental variable 
    @param timeInterval The desired time interval within the environmental variable (i.e. 0 if it is a yearly variable
    or the month index - 0=Jan, 1=Feb etc. - for monthly variables) 
    @return The total of the variable over the whole grid
     */
    double GetEnviroGridTotal(string enviroVariable, unsigned timeInterval) {
        vector<vector< double> > enviroGrid = GetEnviroGrid(enviroVariable, timeInterval);
        double enviroTotal = 0.0;

        for (int ii = 0; ii < NumLatCells; ii += GridCellRarefaction) {
            for (int jj = 0; jj < NumLonCells; jj += GridCellRarefaction) {
                enviroTotal += enviroGrid[ii][ jj];
            }
        }

        return enviroTotal;
    }
    //----------------------------------------------------------------------------------------------
    /** \brief  Return the sum of an environmental variable over specific cells
    @param enviroVariable The environmental variable 
    @param timeInterval The desired time interval within the environmental variable (i.e. 0 if it is a yearly variable
    or the month index - 0=Jan, 1=Feb etc. - for monthly variables) 
    @param cellIndices List of active cells in the model grid 
    @return The total of the variable over the whole grid
     */
    double GetEnviroGridTotal(string enviroVariable, unsigned timeInterval, vector<vector<unsigned>> cellIndices) {
        vector<vector< double> > enviroGrid = GetEnviroGrid(enviroVariable, timeInterval, cellIndices);
        double enviroTotal = 0.0;

        for (int ii = 0; ii < cellIndices.size(); ii++) {
            enviroTotal += enviroGrid[cellIndices[ii][0]][ cellIndices[ii][1]];
        }

        return enviroTotal;
    }
    //----------------------------------------------------------------------------------------------
/** \brief  Get the longitudinal and latitudinal indices of the cell of a viable cell to move to
    @param  latCell The latitudinal index of the focal grid cell 
    @param  lonCell The longitudinal index of the focal grid cell
    @param  v latitudinal displacement
    @param  u longitudinal displacement
    @return The longitudinal and latitudinal cell indices of the cell that lies to the northwest of the focal grid cell
    @remark Currently assumes wrapping in longitude
     */
    const vector<unsigned> getNewCell(const unsigned& latCell, const unsigned& lonCell, const int& v, const int& u) {
        vector<unsigned> Cell = {9999999, 9999999};

        if (latCell + v >= 0 && latCell + v < NumLatCells) {
            int lnc = lonCell + u;
            if (lnc < 0)lnc += NumLonCells;
            if (lnc >= NumLonCells)lnc -= NumLonCells;
                Cell[0] = latCell + v;
                Cell[1] = lnc;
                
        }
        return Cell;
    }
    //----------------------------------------------------------------------------------------------

};
#endif