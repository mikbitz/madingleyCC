#ifndef MODELGRID_H
#define MODELGRID_H
#include <string>
using namespace std;
#include <GridCellCohortHandler.h>
#include <map>
//namespace Madingley
//{
/** \brief
//    /// A class containing the model grid (composed of individual grid cells) along with grid attributes.
//    /// The model grid is referenced by [Lat index, Lon index]\
//    /// <todoD>Check Set and Get state variable methods</todoD>
*/
class ModelGrid
    {
    public:
//        // Private variable to make sure that not more than one grid is instantiated
//        private unsigned NumGrids = 0;
//
//        // Model grid standardised missing value (applied to all grid cells)
//        private double _GlobalMissingValue = -9999;
//
/** \brief
Get the global missing value
*/
//        public double GlobalMissingValue
//        {
//            get { return _GlobalMissingValue; }
//        }
//
//
//        // Field to hold minimum latitude of the grid
//        private float _MinLatitude;    
/** \brief
Get the lower latitude of the lowest cell of the grid
*/
//        public float MinLatitude
//        {
//            get { return _MinLatitude; }
//        }
//
//        // Field to hold minumum longitude of the grid
//        private float _MinLongitude;
/** \brief
Get the leftmost longitude of the leftmost cell of the grid
*/
//        public float MinLongitude
//        {
//            get { return _MinLongitude; }
//        }
//        
//        // Field to hold maximum latitude of the grid
//        private float _MaxLatitude;
/** \brief
Get the lowest latitude of the highest cell in the grid
*/
//        public float MaxLatitude
//        {
//            get { return _MaxLatitude; }
//        }
//        
//        // Field to hold maximum longitude of the grid
//        private float _MaxLongitude;
/** \brief
Get the leftmost longitude of the rightmost cell in the grid
*/
//        public float MaxLongitude
//        {
//            get { return _MaxLongitude; }
//        }
//        
//        // Field to hold latitude resolution of each grid cell
//        private float _LatCellSize;
/** \brief
Get the latitudinal length of each grid cell. Currently assumes all cells are equal sized.
*/
//        public float LatCellSize
//        {
//            get { return _LatCellSize; }
//        }
//        
//        // Field to hold longitude resolution of each grid cell
//        private float _LonCellSize;
/** \brief
Get the longitudinal length of each grid cell. Currently assumes all cells are equal sized. 
*/
//        public float LonCellSize
//        {
//            get { return _LonCellSize; }
//        }
//
/** \brief
The rarefaction of grid cells to be applied to active cells in the model grid
*/
//        private int _GridCellRarefaction;
/** \brief
Get the rarefaction of grid cells to be applied to active cells in the model grid
*/
//        public int GridCellRarefaction
//        { get { return _GridCellRarefaction; } }
//        
/** \brief
The number of latitudinal cells in the model grid
*/
//        private UInt32 _NumLatCells;
/** \brief
Get the number of latitudinal cells in the model grid
*/
//        public UInt32 NumLatCells
//        {
//            get { return _NumLatCells; }
//        }
//        
/** \brief
The number of longitudinal cells in the model grid
*/
//        private UInt32 _NumLonCells;
/** \brief
Get the number of longitudinal cells in the model grid
*/
//        public UInt32 NumLonCells
//        {
//            get { return _NumLonCells; }
//        }
//        
/** \brief
The bottom (southern-most) latitude of each row of grid cells
*/
//        private float[] _Lats;
/** \brief
Get the bottom (southern-most) latitude of each row of grid cells
*/
//        public float[] Lats
//        {
//            get { return _Lats; }
//        }
//        
/** \brief
The left (western-most) longitude of each column of grid cells
*/
//        private float[] _Lons;
/** \brief
Get the left (western-most) longitude of each column of grid cells
*/
//        public float[] Lons
//        {
//            get { return _Lons; }
//        }
//
/** \brief
Array of grid cells
*/
//        GridCell[,] InternalGrid;
//
/** \brief
An array of lists of the functional group indices of each cohort to disperse. Array corresponds to grid cells. The lists correspond to individual cohorts to disperse.
*/
vector< vector <vector<unsigned> > > DeltaFunctionalGroupDispersalArray;
/** \brief
An array of lists of the positions within functional groups of each cohort to disperse. Array corresponds 
to grid cells. The lists correspond to individual cohorts to disperse.
*/
//MB seems clunky - maybe this could be done better! need to get initialisation right...
vector< vector <vector<unsigned> > > DeltaCohortNumberDispersalArray;
//
/** \brief
An array of lists of paired longitude and latitude indices for the grid cells that each cohort will 
to. Array corresponds to grid cells. The lists correspond to paired latitude and longitude indices that 
each cohort will disperse to.
*/
//MB urg!
vector< vector< vector <vector<unsigned> > > > DeltaCellToDisperseToArray;
//        
/** \brief
An array of lists of cells that cohorts in a given grid cell can potentially disperse to (i.e. adjacent cells
in the same realm). Array corresponds to focal grid cells. Lists correspond to cells that cohorts could
disperse to from these focal cells.
*/
vector< vector <vector<unsigned> > > CellsForDispersal;
//
/** \brief
Analagous to the array of lists CellsForDispersal, but instead of containing the identities of the cells that are dispersable to,
instead each array element contains a unsigned list which is coded to correspond to directions:
1. N, 2. NE, 3. E, 4. SE, 5. S, 6. SW, 7 W, 8, NW.
Each item in the list corresponds to the analagous item in CellsForDispersal, and indicates to which direction the cell for dispersal lies. 
This is used by the advective dispersal class, in order to check whether advective dispersal in a particular direction can actually occur.
*/
vector< vector <vector<unsigned> > > CellsForDispersalDirection;
//
/** \brief
The heights of grid cells in each latitudinal band
*/
        vector<double> CellHeightsKm;
/** \brief
The widths of grid cells in each latitudinal band
*/
        vector<double> CellWidthsKm;
//
/** \brief
An instance of the simple random number generator class
*/
//        private NonStaticRNG RandomNumberGenerator = new NonStaticRNG();
//
/** \brief
Instance of the class to perform general functions
*/
//        private UtilityFunctions Utilities;
//
/** \brief
Thread-local variables for tracking extinction and production of cohorts

<todo>Needs a little tidying and checking of access levels</todo>
*/
//        private class ThreadLockedParallelVariablesModelGrid
//        {
//            /// <summary>
//            /// Thread-locked variable to track the cohort ID to assign to newly produced cohorts
//            /// </summary>
//            public Int64 NextCohortIDThreadLocked;
//
//        }
//
/** \brief
Constructor for model grid: assigns grid properties and initialises the grid cells

@param minLat Minimum grid latitude (degrees) 
@param minLon Minimum grid longitude (degrees, currently -180 to 180) 
@param maxLat Maximum grid latitude (degrees) 
@param maxLon Maximum grid longitude (degrees, currently -180 to 180) 
@param latCellSize Latitudinal resolution of grid cell 
@param lonCellSize Longitudinal resolution of grid cell 
@param cellRarefaction The rarefaction to be applied to active grid cells in the model 
@param enviroStack Environmental data layers 
@param cohortFunctionalGroups The functional group definitions for cohorts in the model 
@param stockFunctionalGroups The functional group definitions for stocks in the model 
@param globalDiagnostics Global daignostic variables 
@param nextCohortID The unique ID number to be applied to the next cohort created 
@param tracking Whether process-tracking is enabled 
@param DrawRandomly Whether the model is set to use a random draw 
*/
//        public ModelGrid(float minLat, float minLon,float maxLat,float maxLon,float latCellSize,float lonCellSize, int cellRarefaction, 
//            SortedList<string,EnviroData> enviroStack, FunctionalGroupDefinitions cohortFunctionalGroups, FunctionalGroupDefinitions
//            stockFunctionalGroups, SortedList<string, double> globalDiagnostics, Boolean tracking, Boolean DrawRandomly, Boolean specificLocations)
//        {
//            // Add one to the counter of the number of grids. If there is more than one model grid, exit the program with a debug crash.
//            NumGrids = NumGrids + 1;
//            //Debug.Assert(NumGrids < 2, "You have initialised more than one grid on which to apply models. At present, this is not supported");
//
//            // Initialise the utility functions
//            Utilities = new UtilityFunctions();
//
//            // Seed the random number generator
//            // Set the seed for the random number generator
//            RandomNumberGenerator = new NonStaticRNG();
//            if (DrawRandomly)
//            {
//                RandomNumberGenerator.SetSeedFromSystemTime();
//            }
//            else
//            {
//                RandomNumberGenerator.SetSeed(4315);
//            }
//
//            // CURRENTLY DEFINING MODEL CELLS BY BOTTOM LEFT CORNER
//            _MinLatitude = minLat;
//            _MinLongitude = minLon;
//            _MaxLatitude = maxLat;
//            _MaxLongitude = maxLon;
//            _LatCellSize = latCellSize;
//            _LonCellSize = lonCellSize;
//            _GridCellRarefaction = cellRarefaction;
//
//            // Check to see if the number of grid cells is an integer
//            Debug.Assert((((_MaxLatitude - _MinLatitude) % _LatCellSize) == 0), "Error: number of grid cells is non-integer: check cell size");
//
//            
//            _NumLatCells = (UInt32)((_MaxLatitude - _MinLatitude) / _LatCellSize);
//            _NumLonCells = (UInt32)((_MaxLongitude - _MinLongitude) / _LonCellSize);
//            _Lats = new float[_NumLatCells];
//            _Lons = new float[_NumLonCells];
//
//            // Set up latitude and longitude vectors - lower left
//            for (int ii = 0; ii < _NumLatCells; ii++)
//            {
//                _Lats[ii] = _MinLatitude + ii * _LatCellSize;
//            }
//            for (int jj = 0; jj < _NumLonCells; jj++)
//            {
//                _Lons[jj] = _MinLongitude + jj * _LonCellSize;
//            }
//            
//
//            // Instantiate a grid of grid cells
//            InternalGrid = new GridCell[_NumLatCells, _NumLonCells];
//
//            // Instantiate the arrays of lists of cohorts to disperse
//            DeltaFunctionalGroupDispersalArray = new List<unsigned>[_NumLatCells, _NumLonCells];
//            DeltaCohortNumberDispersalArray = new List<unsigned>[_NumLatCells, _NumLonCells];
//
//            // Instantiate the array of lists of grid cells to disperse those cohorts to
//            DeltaCellToDisperseToArray = new List<unsigned[]>[_NumLatCells, _NumLonCells];
//
//            // An array of lists of cells to which organisms in each cell can disperse to; includes all cells which contribute to the 
//            // perimeter list, plus diagonal cells if they are in the same realm
//            CellsForDispersal = new List<unsigned[]>[_NumLatCells, _NumLonCells];
//
//            // An array of lists of directions corresponding to cells which organisms can disperse to
//            CellsForDispersalDirection = new List<unsigned>[_NumLatCells, _NumLonCells];
//
//            Console.WriteLine("Initialising grid cell environment:");
//
//
//            // Loop through to set up model grid
//            for (int ii = 0; ii < _NumLatCells; ii+=GridCellRarefaction)
//            {
//                for (int jj = 0; jj < _NumLonCells; jj+=GridCellRarefaction)
//                {
//                    InternalGrid[ii, jj] = new GridCell(_Lats[ii],(unsigned)ii, _Lons[jj],(unsigned)jj, LatCellSize, LonCellSize, enviroStack,
//                        GlobalMissingValue, cohortFunctionalGroups, stockFunctionalGroups, globalDiagnostics, tracking, specificLocations);
//                    CellsForDispersal[ii,jj] = new List<unsigned[]>();
//                    CellsForDispersalDirection[ii, jj] = new List<unsigned>();
//                    Console.Write("\rRow {0} of {1}", ii+1, NumLatCells/GridCellRarefaction);
//                }
//            }
//            Console.WriteLine("");
//            Console.WriteLine("");
//
//
//            InterpolateMissingValues();
//
//
//            // Fill in the array of dispersable perimeter lengths for each grid cell
//            CalculatePerimeterLengthsAndCellsDispersableTo();
//
//            CellHeightsKm = new double[_Lats.Length];
//            CellWidthsKm = new double[_Lats.Length];
//
//            // Calculate the lengths of widths of grid cells in each latitudinal strip
//            // Assume that we are at the midpoint of each cell when calculating lengths
//            for (int ii = 0; ii < _Lats.Length; ii++)
//            {
//                 CellHeightsKm[ii] = Utilities.CalculateLengthOfDegreeLatitude(_Lats[ii] + _LatCellSize / 2) * _LatCellSize;
//                 CellWidthsKm[ii] = Utilities.CalculateLengthOfDegreeLongitude(_Lats[ii] + _LatCellSize / 2) * _LonCellSize;
//            }
//        }
//
/** \brief
Overloaded constructor for model grid to construct the grid for specific locations

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
//        public ModelGrid(float minLat, float minLon, float maxLat, float maxLon, float latCellSize, float lonCellSize, List<unsigned[]> cellList, 
//            SortedList<string, EnviroData> enviroStack, FunctionalGroupDefinitions cohortFunctionalGroups,
//            FunctionalGroupDefinitions stockFunctionalGroups, SortedList<string, double> globalDiagnostics, Boolean tracking, 
//            Boolean specificLocations, Boolean runInParallel)
//        { 
//            // Add one to the counter of the number of grids. If there is more than one model grid, exit the program with a debug crash.
//            NumGrids = NumGrids + 1;
//            //Debug.Assert(NumGrids < 2, "You have initialised more than one grid on which to apply models. At present, this is not supported");
//
//            // Initialise the utility functions
//            Utilities = new UtilityFunctions();
//
//            // CURRENTLY DEFINING MODEL CELLS BY BOTTOM LEFT CORNER
//            _MinLatitude = minLat;
//            _MinLongitude = minLon;
//            _MaxLatitude = maxLat;
//            _MaxLongitude = maxLon;
//            _LatCellSize = latCellSize;
//            _LonCellSize = lonCellSize;
//            _GridCellRarefaction = 1;
//
//            // Check to see if the number of grid cells is an integer
//            Debug.Assert((((_MaxLatitude - _MinLatitude) % _LatCellSize) == 0), "Error: number of grid cells is non-integer: check cell size");
//
//
//            _NumLatCells = (UInt32)((_MaxLatitude - _MinLatitude) / _LatCellSize);
//            _NumLonCells = (UInt32)((_MaxLongitude - _MinLongitude) / _LonCellSize);
//            _Lats = new float[_NumLatCells];
//            _Lons = new float[_NumLonCells];
//
//            // Set up latitude and longitude vectors - lower left
//            for (int ii = 0; ii < _NumLatCells; ii++)
//            {
//                _Lats[ii] = _MinLatitude + ii * _LatCellSize;
//            }
//            for (int jj = 0; jj < _NumLonCells; jj++)
//            {
//                _Lons[jj] = _MinLongitude + jj * _LonCellSize;
//            }
//
//
//            // Set up a grid of grid cells
//            InternalGrid = new GridCell[_NumLatCells, _NumLonCells];
//
//            // Instantiate the arrays of lists of cohorts to disperse
//            DeltaFunctionalGroupDispersalArray = new List<unsigned>[_NumLatCells, _NumLonCells];
//            DeltaCohortNumberDispersalArray = new List<unsigned>[_NumLatCells, _NumLonCells];
//
//            // Instantiate the array of lists of grid cells to disperse those cohorts to
//            DeltaCellToDisperseToArray = new List<unsigned[]>[_NumLatCells, _NumLonCells];
//
//
//            // An array of lists of cells to which organisms in each cell can disperse to; includes all cells which contribute to the 
//            // perimeter list, plus diagonal cells if they are in the same realm
//            CellsForDispersal = new List<unsigned[]>[_NumLatCells, _NumLonCells];
//
//            // An array of lists of directions corresponding to cells which organisms can disperse to
//            CellsForDispersalDirection = new List<unsigned>[_NumLatCells, _NumLonCells];
//
//            Console.WriteLine("Initialising grid cell environment:");
//
//
//            int Count = 0;
//
//            int NCells = cellList.Count;
//
//            if (!runInParallel)
//            {
//                // Loop over cells to set up the model grid
//                for (int ii = 0; ii < cellList.Count; ii++)
//                {
//                    // Create the grid cell at the specified position
//                    InternalGrid[cellList[ii][0], cellList[ii][1]] = new GridCell(_Lats[cellList[ii][0]], cellList[ii][0],
//                        _Lons[cellList[ii][1]], cellList[ii][1], latCellSize, lonCellSize, enviroStack, _GlobalMissingValue,
//                        cohortFunctionalGroups, stockFunctionalGroups, globalDiagnostics, tracking, specificLocations);
//                    if (!specificLocations)
//                    {
//                        CellsForDispersal[cellList[ii][0], cellList[ii][1]] = new List<unsigned[]>();
//                        CellsForDispersalDirection[cellList[ii][0], cellList[ii][1]] = new List<unsigned>();
//                        Count++;
//                        Console.Write("\rInitialised {0} of {1}", Count, NCells);
//                    }
//                    else
//                    {
//                        Console.Write("\rRow {0} of {1}", ii + 1, NumLatCells / GridCellRarefaction);
//                        Console.WriteLine("");
//                        Console.WriteLine("");
//                    }
//                }
//            }
//            else
//            {
//
//                // Run a parallel loop over rows
//                Parallel.For(0, NCells, ii =>
//                {
//                    // Create the grid cell at the specified position
//                    InternalGrid[cellList[ii][0], cellList[ii][1]] = new GridCell(_Lats[cellList[ii][0]], cellList[ii][0],
//                        _Lons[cellList[ii][1]], cellList[ii][1], latCellSize, lonCellSize, enviroStack, _GlobalMissingValue,
//                        cohortFunctionalGroups, stockFunctionalGroups, globalDiagnostics, tracking, specificLocations);
//                    if (!specificLocations)
//                    {
//                        CellsForDispersal[cellList[ii][0], cellList[ii][1]] = new List<unsigned[]>();
//                        CellsForDispersalDirection[cellList[ii][0], cellList[ii][1]] = new List<unsigned>();
//                    }
//
//                    Count++;
//                    Console.Write("\rInitialised {0} of {1}", Count, NCells);
//                }
//                 );
//
//            }
//
//
//            if (!specificLocations)
//            {
//                InterpolateMissingValues();
//
//
//                // Fill in the array of dispersable perimeter lengths for each grid cell
//                CalculatePerimeterLengthsAndCellsDispersableTo();
//
//                CellHeightsKm = new double[_Lats.Length];
//                CellWidthsKm = new double[_Lats.Length];
//
//                // Calculate the lengths of widths of grid cells in each latitudinal strip
//                // Assume that we are at the midpoint of each cell when calculating lengths
//                for (int ii = 0; ii < _Lats.Length; ii++)
//                {
//                    CellHeightsKm[ii] = Utilities.CalculateLengthOfDegreeLatitude(_Lats[ii] + _LatCellSize / 2) * _LatCellSize;
//                    CellWidthsKm[ii] = Utilities.CalculateLengthOfDegreeLongitude(_Lats[ii] + _LatCellSize / 2) * _LonCellSize;
//                }
//            }
//
//            Console.WriteLine("\n");
//
//        }
//
/** \brief
Estimates missing environmental data for grid cells by interpolation
*/
//        public void InterpolateMissingValues()
//        {
//            SortedList<string, double[]> WorkingCellEnvironment = new SortedList<string, double[]>();
//            Boolean Changed = false;
//
//            for (unsigned ii = 0; ii < _NumLatCells; ii++)
//            {
//                for (unsigned jj = 0; jj < _NumLonCells; jj++)
//                {
//                    WorkingCellEnvironment = GetCellEnvironment(ii, jj);
//
//                    // If the cell environment does not contain valid NPP data then interpolate values
//                    if (!InternalGrid[ii, jj].ContainsData(WorkingCellEnvironment["NPP"], WorkingCellEnvironment["Missing Value"][0]))
//                    {
//                        //If NPP doesn't exist the interpolate from surrounding values (of the same realm)
//                        WorkingCellEnvironment["NPP"] = GetInterpolatedValues(ii, jj, GetCellLatitude(ii), GetCellLongitude(jj), "NPP", WorkingCellEnvironment["Realm"][0]);
//                        
//                        //Calculate NPP seasonality - for use in converting annual NPP estimates to monthly
//                        WorkingCellEnvironment["Seasonality"] = InternalGrid[ii, jj].CalculateNPPSeasonality(WorkingCellEnvironment["NPP"], WorkingCellEnvironment["Missing Value"][0]);
//                        Changed = true;
//                    }
//                    // Otherwise convert the missing data values to zeroes where they exist amongst valid data eg in polar regions.
//                    else
//                    {
//                        WorkingCellEnvironment["NPP"] = InternalGrid[ii, jj].ConvertMissingValuesToZero(WorkingCellEnvironment["NPP"], WorkingCellEnvironment["Missing Value"][0]);
//                    }
//
//                    // If the cell environment does not contain valid monthly mean diurnal temperature range data then interpolate values
//                    if (InternalGrid[ii, jj].ContainsMissingValue(WorkingCellEnvironment["DiurnalTemperatureRange"], WorkingCellEnvironment["Missing Value"][0]))
//                    {
//                        //If NPP doesn't exist the interpolate from surrounding values (of the same realm)
//                        WorkingCellEnvironment["DiurnalTemperatureRange"] = FillWithInterpolatedValues(ii, jj, GetCellLatitude(ii), GetCellLongitude(jj), "DiurnalTemperatureRange", WorkingCellEnvironment["Realm"][0]);
//
//                        Changed = true;
//                    }
//
//                    // Same for u and v velocities
//                    if (!InternalGrid[ii, jj].ContainsData(WorkingCellEnvironment["uVel"], WorkingCellEnvironment["Missing Value"][0]))
//                    {
//                        //If u doesn't exist the interpolate from surrounding values (of the same realm)
//                        WorkingCellEnvironment["uVel"] = GetInterpolatedValues(ii, jj, GetCellLatitude(ii), GetCellLongitude(jj), "uVel", WorkingCellEnvironment["Realm"][0]);
//
//                        Changed = true;
//                    }
//                    // Otherwise convert the missing data values to zeroes where they exist amongst valid data eg in polar regions.
//                    else
//                    {
//                        WorkingCellEnvironment["uVel"] = InternalGrid[ii, jj].ConvertMissingValuesToZero(WorkingCellEnvironment["uVel"], WorkingCellEnvironment["Missing Value"][0]);
//                    }
//
//                    if (!InternalGrid[ii, jj].ContainsData(WorkingCellEnvironment["vVel"], WorkingCellEnvironment["Missing Value"][0]))
//                    {
//                        //If v vel doesn't exist the interpolate from surrounding values (of the same realm)
//                        WorkingCellEnvironment["vVel"] = GetInterpolatedValues(ii, jj, GetCellLatitude(ii), GetCellLongitude(jj), "vVel", WorkingCellEnvironment["Realm"][0]);
//
//                        Changed = true;
//                    }
//                    // Otherwise convert the missing data values to zeroes where they exist amongst valid data eg in polar regions.
//                    else
//                    {
//                        WorkingCellEnvironment["vVel"] = InternalGrid[ii, jj].ConvertMissingValuesToZero(WorkingCellEnvironment["vVel"], WorkingCellEnvironment["Missing Value"][0]);
//                    }
//                    
//                    if(Changed) InternalGrid[ii, jj].CellEnvironment = WorkingCellEnvironment;
//                }
//            }
//        }
//
/** \brief
Calculate the weighted average of surrounding grid cell data, where those grid cells are of the specified realm and contain
non missing data values

@param latIndex Index of the latitude cell for which the weighted average over surrounding cells is requested 
@param lonIndex Index of the longitude cell for which the weighted average over surrounding cells is requested 
@param lat Latitude of the cell for which the weighted value is requested 
@param lon Longitude of the cell for which the weighted value is requested 
@param dataName Names of the data for which weighted value is requested 
@param realm Realm of the grid cell for which data is to be averaged over 
@return The weighted average value of the specified data type across surrounding grid cells of the specified realm
*/
//        private double[] GetInterpolatedValues(unsigned latIndex, unsigned lonIndex, double lat, double lon, string dataName, double realm)
//        {
//            SortedList<string, double[]> TempCellEnvironment = GetCellEnvironment(latIndex, lonIndex);
//            double[] InterpData = new double[TempCellEnvironment[dataName].Length];
//            unsigned[] InterpCount = new unsigned[TempCellEnvironment[dataName].Length];
//
//            unsigned LowerLatIndex = latIndex - 1;
//            unsigned UpperLatIndex = latIndex + 1;
//            unsigned LowerLonIndex = lonIndex - 1;
//            unsigned UpperLonIndex = lonIndex + 1;
//
//
//            if (latIndex == 0) LowerLatIndex = latIndex;
//            if (lat.CompareTo(this.MaxLatitude) == 0) UpperLatIndex = latIndex;
//
//            if (lonIndex == 0) LowerLonIndex = lonIndex;
//            if (lon.CompareTo(this.MaxLongitude) == 0) UpperLonIndex = lonIndex;
//
//            //Loop over surrounding cells in the datalayer
//            for (unsigned ii = LowerLatIndex; ii <= UpperLatIndex; ii++)
//            {
//                for (unsigned jj = LowerLonIndex; jj < UpperLonIndex; jj++)
//                {
//                    if (ii < _NumLatCells && jj < _NumLonCells)
//                    {
//                        TempCellEnvironment = GetCellEnvironment(ii, jj);
//
//                        for (unsigned hh = 0; hh < InterpData.Length; hh++)
//                        {
//                            //If the cell contains data then sum this and increment count
//                            if (TempCellEnvironment[dataName][hh] != TempCellEnvironment["Missing Value"][0] && TempCellEnvironment["Realm"][0] == realm)
//                            {
//                                InterpData[hh] += TempCellEnvironment[dataName][hh];
//                                InterpCount[hh]++;
//                            }
//                        }
//                    }
//                }
//            }
//
//            //take the mean over surrounding valid cells for each timestep
//            for (int hh = 0; hh < InterpData.Length; hh++)
//            {
//                if (InterpCount[hh] > 0)
//                {
//                    InterpData[hh] /= InterpCount[hh];
//                }
//                else
//                {
//                    InterpData[hh] = 0.0;
//                }
//            }
//            return InterpData;
//        }
//
/** \brief
Calculate the weighted average of surrounding grid cell data, where those grid cells are of the specified realm and contain
non missing data values

@param latIndex Index of the latitude cell for which the weighted average over surrounding cells is requested 
@param lonIndex Index of the longitude cell for which the weighted average over surrounding cells is requested 
@param lat Latitude of the cell for which the weighted value is requested 
@param lon Longitude of the cell for which the weighted value is requested 
@param dataName Names of the data for which weighted value is requested 
@param realm Realm of the grid cell for which data is to be averaged over 
@returns The weighted average value of the specified data type across surrounding grid cells of the specified realm
*/
//        private double[] FillWithInterpolatedValues(unsigned latIndex, unsigned lonIndex, double lat, double lon, string dataName, double realm)
//        {
//            SortedList<string, double[]> TempCellEnvironment = GetCellEnvironment(latIndex, lonIndex);
//            double[] InterpData = new double[TempCellEnvironment[dataName].Length];
//            unsigned[] InterpCount = new unsigned[TempCellEnvironment[dataName].Length];
//            unsigned LowerLatIndex = latIndex - 1;
//            unsigned UpperLatIndex = latIndex + 1;
//            unsigned LowerLonIndex = lonIndex - 1;
//            unsigned UpperLonIndex = lonIndex + 1;
//
//
//            if (latIndex == 0) LowerLatIndex = latIndex;
//            if (lat.CompareTo(this.MaxLatitude) == 0) UpperLatIndex = latIndex;
//
//            if (lonIndex == 0) LowerLonIndex = lonIndex;
//            if (lon.CompareTo(this.MaxLongitude) == 0) UpperLonIndex = lonIndex;
//
//            for (unsigned hh = 0; hh < InterpData.Length; hh++)
//            {
//                if (TempCellEnvironment[dataName][hh] == TempCellEnvironment["Missing Value"][0])
//                {
//                    //Loop over surrounding cells in the datalayer
//                    for (unsigned ii = LowerLatIndex; ii <= UpperLatIndex; ii++)
//                    {
//                        for (unsigned jj = LowerLonIndex; jj <= UpperLonIndex; jj++)
//                        {
//                            if (ii < _NumLatCells && jj < _NumLonCells)
//                            {
//                                TempCellEnvironment = GetCellEnvironment(ii, jj);
//
//                                //If the cell contains data then sum this and increment count
//                                if (TempCellEnvironment[dataName][hh] != TempCellEnvironment["Missing Value"][0] && TempCellEnvironment["Realm"][0] == realm)
//                                {
//                                    InterpData[hh] += TempCellEnvironment[dataName][hh];
//                                    InterpCount[hh]++;
//                                }
//
//                            }
//                        }
//                    }
//                    //take the mean over surrounding valid cells for each timestep
//                    if (InterpCount[hh] > 0)
//                    {
//                        InterpData[hh] /= InterpCount[hh];
//                    }
//                    else
//                    {
//                        InterpData[hh] = 0.0;
//                    }
//                }
//                else
//                {
//                    InterpData[hh] = TempCellEnvironment[dataName][hh];
//                }
//            }
//
//
//            return InterpData;
//        }
//
/** \brief
Seed the stocks and cohorts for all active cells in the model grid

@param cellIndices A list of the active cells in the model grid 
@param cohortFunctionalGroupDefinitions The functional group definitions for cohorts in the model 
@param stockFunctionalGroupDefinitions The functional group definitions for stocks in the model 
@param globalDiagnostics A list of global diagnostic variables 
@param nextCohortID The ID number to be assigned to the next produced cohort 
@param tracking Whether process-tracking is enabled 
@param DrawRandomly Whether the model is set to use a random draw 
@param dispersalOnly Whether to run dispersal only (i.e. to turn off all other ecological processes 
@param processTrackers An instance of the ecological process tracker 
*/
//        public void SeedGridCellStocksAndCohorts(List<unsigned[]> cellIndices, FunctionalGroupDefinitions cohortFunctionalGroupDefinitions, 
//            FunctionalGroupDefinitions stockFunctionalGroupDefinitions, SortedList<string,double> globalDiagnostics, ref Int64 nextCohortID,
//            Boolean tracking, Boolean DrawRandomly, Boolean dispersalOnly, string dispersalOnlyType)
//        {
//            int ii = 1;
//            Console.WriteLine("Seeding grid cell stocks and cohorts:");
//
//            //Work out how many cohorts are to be seeded in each grid cell - split by realm as different set of cohorts initialised by realm
//            double TotalTerrestrialCellCohorts = 0;
//            double TotalMarineCellCohorts = 0;
//
//            int[] TerrestrialFunctionalGroups = cohortFunctionalGroupDefinitions.GetFunctionalGroupIndex("Realm", "Terrestrial", false);
//            if (TerrestrialFunctionalGroups == null)
//            {
//                TotalTerrestrialCellCohorts = 0;
//            }
//            else
//            {
//                foreach (int F in TerrestrialFunctionalGroups)
//                {
//                    TotalTerrestrialCellCohorts += cohortFunctionalGroupDefinitions.GetBiologicalPropertyOneFunctionalGroup("Initial number of GridCellCohorts", F);
//                }
//            }
//
//
//            int[] MarineFunctionalGroups = cohortFunctionalGroupDefinitions.GetFunctionalGroupIndex("Realm", "Marine", false);
//            if (MarineFunctionalGroups == null)
//            {
//                TotalMarineCellCohorts = 0;
//            }
//            else
//            {
//                foreach (int F in MarineFunctionalGroups)
//                {
//                    TotalMarineCellCohorts += cohortFunctionalGroupDefinitions.GetBiologicalPropertyOneFunctionalGroup("Initial number of GridCellCohorts", F);
//                }
//            }
//
//            foreach (unsigned[] cellIndexPair in cellIndices)
//            {
//                if (dispersalOnly)
//                {
//                    if (dispersalOnlyType == "diffusion")
//                    {
//                        // Diffusive dispersal
//
//                        if ((cellIndexPair[0] == 90) && (cellIndexPair[1] == 180))
//                        {
//                            InternalGrid[cellIndexPair[0], cellIndexPair[1]].SeedGridCellCohortsAndStocks(cohortFunctionalGroupDefinitions,
//                            stockFunctionalGroupDefinitions, globalDiagnostics, ref nextCohortID, tracking, TotalTerrestrialCellCohorts, TotalMarineCellCohorts,
//                            DrawRandomly, false);
//                        }
//                        else if ((cellIndexPair[0] == 95) && (cellIndexPair[1] == 110))
//                        {
//                            InternalGrid[cellIndexPair[0], cellIndexPair[1]].SeedGridCellCohortsAndStocks(cohortFunctionalGroupDefinitions,
//                            stockFunctionalGroupDefinitions, globalDiagnostics, ref nextCohortID, tracking, TotalTerrestrialCellCohorts, TotalMarineCellCohorts,
//                            DrawRandomly, false);
//                        }
//                        else
//                        {
//                            InternalGrid[cellIndexPair[0], cellIndexPair[1]].SeedGridCellCohortsAndStocks(cohortFunctionalGroupDefinitions,
//                            stockFunctionalGroupDefinitions, globalDiagnostics, ref nextCohortID, tracking, TotalTerrestrialCellCohorts, TotalMarineCellCohorts,
//                            DrawRandomly, true);
//                        }
//                        Console.Write("\rGrid Cell: {0} of {1}", ii++, cellIndices.Count);
//                    }
//                    else if (dispersalOnlyType == "advection")
//                    {
//                        // Advective dispersal
//                        /*
//                        if ((cellIndexPair[0] == 58) && (cellIndexPair[1] == 225))
//                        {
//                            InternalGrid[cellIndexPair[0], cellIndexPair[1]].SeedGridCellCohortsAndStocks(cohortFunctionalGroupDefinitions,
//                            stockFunctionalGroupDefinitions, globalDiagnostics, ref nextCohortID, tracking, TotalTerrestrialCellCohorts, TotalMarineCellCohorts,
//                            DrawRandomly, false);
//                        }
//                        else if ((cellIndexPair[0] == 95) && (cellIndexPair[1] == 110))
//                        {
//                            InternalGrid[cellIndexPair[0], cellIndexPair[1]].SeedGridCellCohortsAndStocks(cohortFunctionalGroupDefinitions,
//                            stockFunctionalGroupDefinitions, globalDiagnostics, ref nextCohortID, tracking, TotalTerrestrialCellCohorts, TotalMarineCellCohorts,
//                            DrawRandomly, false);
//                        }
//                        else
//                        {
//                            InternalGrid[cellIndexPair[0], cellIndexPair[1]].SeedGridCellCohortsAndStocks(cohortFunctionalGroupDefinitions,
//                            stockFunctionalGroupDefinitions, globalDiagnostics, ref nextCohortID, tracking, TotalTerrestrialCellCohorts, TotalMarineCellCohorts,
//                            DrawRandomly, true);
//                        }
//                        */
//                        if (InternalGrid[cellIndexPair[0], cellIndexPair[1]].CellEnvironment["Realm"][0] == 1.0)
//                        {
//                            InternalGrid[cellIndexPair[0], cellIndexPair[1]].SeedGridCellCohortsAndStocks(cohortFunctionalGroupDefinitions,
//    stockFunctionalGroupDefinitions, globalDiagnostics, ref nextCohortID, tracking, TotalTerrestrialCellCohorts, TotalMarineCellCohorts,
//    DrawRandomly, true);
//                        }
//                        else
//                        {
//                            InternalGrid[cellIndexPair[0], cellIndexPair[1]].SeedGridCellCohortsAndStocks(cohortFunctionalGroupDefinitions,
//    stockFunctionalGroupDefinitions, globalDiagnostics, ref nextCohortID, tracking, TotalTerrestrialCellCohorts, TotalMarineCellCohorts,
//    DrawRandomly, false);
//                        }
//                        Console.Write("\rGrid Cell: {0} of {1}", ii++, cellIndices.Count);
//                    }
//                    else if (dispersalOnlyType == "responsive")
//                    {
//                        // Responsive dispersal
//
//                        InternalGrid[cellIndexPair[0], cellIndexPair[1]].SeedGridCellCohortsAndStocks(cohortFunctionalGroupDefinitions,
//                        stockFunctionalGroupDefinitions, globalDiagnostics, ref nextCohortID, tracking, TotalTerrestrialCellCohorts, TotalMarineCellCohorts,
//                        DrawRandomly, true);
//
//                        Console.Write("\rGrid Cell: {0} of {1}", ii++, cellIndices.Count);
//                    }
//                    else
//                    {
//                        Debug.Fail("Dispersal only type not recognized from initialisation file");
//                    }
//                }
//
//                else
//                {
//                    InternalGrid[cellIndexPair[0], cellIndexPair[1]].SeedGridCellCohortsAndStocks(cohortFunctionalGroupDefinitions,
//                        stockFunctionalGroupDefinitions, globalDiagnostics, ref nextCohortID, tracking, TotalTerrestrialCellCohorts, TotalMarineCellCohorts,
//                        DrawRandomly, false);
//                    Console.Write("\rGrid Cell: {0} of {1}", ii++, cellIndices.Count);
//                }
//            }
//            Console.WriteLine("");
//            Console.WriteLine("");
//        }
//
/** \brief
Returns the stocks within the specified grid cell

@param latIndex Latitude index 
@param lonIndex Longitude index 
@returns The stock handler for the specified grid cell
*/
//        public GridCellStockHandler GetGridCellStocks(unsigned latIndex, unsigned lonIndex)
//        {
//            return InternalGrid[latIndex, lonIndex].GridCellStocks;
//        }
//
/** \brief
Sets the stocks in the specified grid cell to the passed stocks

@param newGridCellStocks New stocks for the grid cell 
@param latIndex Latitude index 
@param lonIndex Longitude index 
*/
//        public void SetGridCellStocks(GridCellStockHandler newGridCellStocks, unsigned latIndex, unsigned lonIndex)
//        {
//            InternalGrid[latIndex, lonIndex].GridCellStocks = newGridCellStocks;
//        }
//
/** \brief
Returns the array (indexed by functional group) of lists of gridCellCohorts for the specified grid cell
@param latIndex Latitude index of grid cell 
@param lonIndex Longitude index of grid cell 
@returns Array (indexed by functional group) of lists of gridCellCohorts
*/
//MB This type of function should return a reference?? Don't want to copy surely?
       GridCellCohortHandler GetGridCellCohorts(unsigned latIndex, unsigned lonIndex)
        {
//            return InternalGrid[latIndex, lonIndex].GridCellCohorts;
        }
//
/** \brief
Extracts an individual cohort from a particular grid cell

@param latIndex Latitude index of grid cell 
@param lonIndex Longitude index of grid cell 
@param functionalGroup Functional group of cohort 
@param positionInList Index of cohort position in the list 
@return what?
*/
//        public Cohort GetGridCellIndividualCohort(unsigned latIndex, unsigned lonIndex, int functionalGroup, int positionInList)
//        {
//            return InternalGrid[latIndex, lonIndex].GridCellCohorts[functionalGroup].ElementAt(positionInList);
//        }
//
//        // NOTE TO SELF: These need more error checking, and also the access levels more tightly controlled
/** \brief Remove an individual cohort from a functionall group; necessary due to dispersal moving cohorts from one cell to another

@param latIndex Grid cell latitude index 
@param lonIndex Grid cell longitude index 
@param functionalGroup Cohort functional group 
@param positionInList Position of cohort in the list of that functional group 
*/
//        public void DeleteGridCellIndividualCohort(unsigned latIndex, unsigned lonIndex, int functionalGroup, int positionInList)
//        {
//            InternalGrid[latIndex, lonIndex].GridCellCohorts[functionalGroup].RemoveAt(positionInList);
//        }
//
/** \brief Delete a specified list of cohorts from a grid cell

@param latIndex The latitudinal index of the grid cell to delete cohorts from 
@param lonIndex The longitudinal index of the grid cell to delete cohorts from 
@param cohortFGsToDelete A list of the functional groups that each cohort to delete belongs to 
@param cohortNumbersToDelete A list of the positions with each functional group that each cohort to delete occupies 
\remarks This is inefficient and needs double-checking for errors
*/
//        public void DeleteGridCellIndividualCohorts(unsigned latIndex, unsigned lonIndex, List<unsigned> cohortFGsToDelete, List<unsigned> cohortNumbersToDelete)
//        {
//            
//            // Get the unique functional groups that have cohorts to be removed
//            unsigned[] TempList = cohortFGsToDelete.Distinct().ToArray();
//
//            // Loop over these unique functional  groups
//            for (int ii = 0; ii < TempList.Length; ii++)
//			{
//                // Get the functional group index of the current functional group
//                int FG = (int)TempList[ii];
//
//                // Create a local list to hold the positions of the cohorts to delete from this functional group
//			    List<unsigned> CohortIndexList = new List<unsigned>();
//                // Loop over all cohorts to be deleted
//                for (int jj = 0; jj < cohortFGsToDelete.Count; jj++)
//			    {
//                    // Check whether the functional group correpsonds with the functional group currently being processed 
//                    if (cohortFGsToDelete.ElementAt((int)jj) == FG)
//                    {
//                        // Add the cohort to the list of cohorts to delete
//                        CohortIndexList.Add(cohortNumbersToDelete[jj]);
//                    }
//			    }
//
//                // Sort the list of positions of the cohorts to delete in this functional group
//                CohortIndexList.Sort();
//                // Reverse the list so that the highest positions come first
//                CohortIndexList.Reverse();
//
//                // Loop over cohorts and delete in turn, starting with cohorts in the highest positions
//                for (int kk = 0; kk < CohortIndexList.Count; kk++)
//			    {
//			        InternalGrid[latIndex, lonIndex].GridCellCohorts[FG].RemoveAt((int)CohortIndexList[kk]);                  
//			    }
//
//	        }
//            
//        }
//
/** \brief Replace the gridCellCohorts in a grid cell with a new list of gridCellCohorts

@param newGridCellCohorts The new list of gridCellCohorts 
@param latIndex Grid cell latitude index 
@param lonIndex Grid cell longitude index 
*/
//        public void SetGridCellCohorts(GridCellCohortHandler newGridCellCohorts, unsigned latIndex, unsigned lonIndex)
//        {
//            InternalGrid[latIndex, lonIndex].GridCellCohorts = newGridCellCohorts;
//        }
//
/** \brief
Add a new cohort to an existing list of cohorts in the grid cell - or create a new list if there is not one present

@param latIndex Latitude index of the grid cell 
@param lonIndex Longitude index of the grid cell 
@param functionalGroup Functional group of the cohort (i.e. array index) 
@param cohortToAdd The cohort object to add 
*/
//        public void AddNewCohortToGridCell(unsigned latIndex, unsigned lonIndex, int functionalGroup, Cohort cohortToAdd)
//        {
//            InternalGrid[latIndex, lonIndex].GridCellCohorts[functionalGroup].Add(cohortToAdd);
//        }
//
/** \brief
Return the value of a specified environmental layer from an individual grid cell

@param variableName The name of the environmental lyaer 
@param timeInterval The desired time interval within the environmental variable (i.e. 0 if it is a yearly variable
or the month index - 0=Jan, 1=Feb etc. - for monthly variables) 
@param latCellIndex The latitudinal cell index 
@param lonCellIndex The longitudinal cell index 
@param variableExists Returns false if the environmental layer does not exist, true if it does 
@returns The value of the environmental layer, or a missing value if the environmental layer does not exist
*/
//MB This type of function should return a reference?? Don't want to copy surely?

double GetEnviroLayer(string variableName, unsigned timeInterval, unsigned latCellIndex, unsigned lonCellIndex, bool variableExists)
        {
//            return InternalGrid[latCellIndex, lonCellIndex].GetEnviroLayer(variableName,timeInterval, out variableExists);
        }
//
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
//        public bool SetEnviroLayer(string variableName, unsigned timeInterval, double setValue, unsigned latCellIndex, unsigned lonCellIndex)
//        {
//            return InternalGrid[latCellIndex, lonCellIndex].SetEnviroLayer(variableName,timeInterval, setValue);
//        }
//
//        
/** \brief
Set the value of a given delta type for the specified ecological process within the specified grid cell

@param deltaType The type of delta value to set (e.g. 'biomass', 'abundance' etc.) 
@param ecologicalProcess The name of the ecological process to set the value of delta for 
@param setValue The value to set 
@param latCellIndex The latitudinal index of the cell 
@param lonCellIndex The longitudinal index of the cell 
@returns True if the value is set successfully, false otherwise
*/
//        public bool SetDeltas(string deltaType, string ecologicalProcess, double setValue, unsigned latCellIndex, unsigned lonCellIndex)
//        {
//            return InternalGrid[latCellIndex, lonCellIndex].SetDelta(deltaType, ecologicalProcess, setValue);
//        }
//        
/** \brief
Get the total of a state variable for specific cells

@param variableName The name of the variable 
@param functionalGroups A vector of functional group indices to consider 
@param cellIndices List of indices of active cells in the model grid 
@param stateVariableType A string indicating the type of state variable; 'cohort' or 'stock' 
@returns Summed value of variable over whole grid
<todo>Overload to work with vector and array state variables</todo>
*/
//        public double StateVariableGridTotal(string variableName, string traitValue, int[] functionalGroups, List<unsigned[]> cellIndices, string stateVariableType, MadingleyModelInitialisation initialisation)
//        {
//
//            double tempVal = 0;
//
//            double[,] TempStateVariable = this.GetStateVariableGrid(variableName, traitValue, functionalGroups, cellIndices, stateVariableType, initialisation);
//
//            // Loop through and sum values across a grid, excluding missing values
//            for (int ii = 0; ii < cellIndices.Count; ii++)
//            {
//                tempVal += TempStateVariable[cellIndices[ii][0], cellIndices[ii][1]];
//            }
//
//            return tempVal;
//        }
//
/** \brief
Gets a state variable for specified functional groups of specified entity types in a specified grid cell

@param variableName The name of the variable to get: 'biomass' or 'abundance' 
@param functionalGroups The functional group indices to get the state variable for 
@param latCellIndex The latitudinal index of the cell 
@param lonCellIndex The longitudinal index of the cell 
@param stateVariableType The type of entity to return the state variable for: 'stock' or 'cohort' 
@returns>The state variable for specified functional groups of specified entity types in a specified grid cell
*/
//        public double GetStateVariable(string variableName, string traitValue, int[] functionalGroups, unsigned latCellIndex, unsigned lonCellIndex, string stateVariableType, MadingleyModelInitialisation modelInitialisation)
//        {
//
//            double returnValue = 0.0;
//
//            switch (stateVariableType.ToLower())
//            {
//                case "cohort":
//                    
//                    GridCellCohortHandler TempCohorts = InternalGrid[latCellIndex, lonCellIndex].GridCellCohorts;
//
//                    switch (variableName.ToLower())
//                    {
//                        case "biomass":
//                            if (traitValue != "Zooplankton")
//                            {
//                                foreach (int f in functionalGroups)
//                                {
//                                    foreach (var item in TempCohorts[f])
//                                    {
//                                        returnValue += ((item.IndividualBodyMass + item.IndividualReproductivePotentialMass) * item.CohortAbundance);
//                                    }
//                                }
//                            }
//                            else
//                            {
//                                foreach (int f in functionalGroups)
//                                {
//                                    foreach (var item in TempCohorts[f])
//                                    {
//                                        if (item.IndividualBodyMass <= modelInitialisation.PlanktonDispersalThreshold)
//                                        returnValue += ((item.IndividualBodyMass + item.IndividualReproductivePotentialMass) * item.CohortAbundance);
//                                    }
//                                }
//                            }
//                            break;
//
//                        case "abundance":
//                            if (traitValue != "Zooplankton")
//                            {
//                                foreach (int f in functionalGroups)
//                                {
//                                    foreach (var item in TempCohorts[f])
//                                    {
//                                        returnValue += item.CohortAbundance;
//                                    }
//                                }
//                            }
//                            else
//                            {
//                                foreach (int f in functionalGroups)
//                                {
//                                    foreach (var item in TempCohorts[f])
//                                    {
//                                        if (item.IndividualBodyMass <= modelInitialisation.PlanktonDispersalThreshold)
//                                        returnValue += item.CohortAbundance;
//                                    }
//                                }
//                            }
//                            break;
//
//                        default:
//                            Debug.Fail("For cohorts, state variable name must be either 'biomass' or 'abundance'");
//                            break;
//                    }
//                    break;
//
//                case "stock":
//                    GridCellStockHandler TempStocks = InternalGrid[latCellIndex, lonCellIndex].GridCellStocks;
//
//                    switch (variableName.ToLower())
//                    {
//                        case "biomass":
//                            foreach (int f in functionalGroups)
//                            {
//                                foreach (var item in TempStocks[f])
//                                {
//                                    returnValue += item.TotalBiomass;
//                                }
//                            }
//                            break;
//                        default:
//                            Debug.Fail("For stocks, state variable name must be 'biomass'");
//                            break;
//                    }
//                    break;
//
//                default:
//                    Debug.Fail("State variable type must be either 'cohort' or 'stock'");
//                    break;
//
//            }
//
//            
//
//
//
//
//
//            return returnValue;
//        }
//
/** \brief
Gets a state variable density for specified functional groups of specified entity types in a specified grid cell

@param variableName The name of the variable to get: 'biomass' or 'abundance' 
@param functionalGroups The functional group indices to get the state variable for 
@param latCellIndex The latitudinal index of the cell 
@param lonCellIndex The longitudinal index of the cell 
@param stateVariableType The type of entity to return the state variable for: 'stock' or 'cohort' 
\returns The state variable density for specified functional groups of specified entity types in a specified grid cell
*/
//        public double GetStateVariableDensity(string variableName, string traitValue, int[] functionalGroups, unsigned latCellIndex, unsigned lonCellIndex, string stateVariableType, MadingleyModelInitialisation modelInitialisation)
//        {
//
//            double returnValue = 0.0;
//
//            switch (stateVariableType.ToLower())
//            {
//                case "cohort":
//
//                    GridCellCohortHandler TempCohorts = InternalGrid[latCellIndex, lonCellIndex].GridCellCohorts;
//
//                    switch (variableName.ToLower())
//                    {
//                        case "biomass":
//                            if (traitValue != "Zooplankton (all)")
//                            {
//                                foreach (int f in functionalGroups)
//                                {
//                                    foreach (var item in TempCohorts[f])
//                                    {
//                                        returnValue += ((item.IndividualBodyMass + item.IndividualReproductivePotentialMass) * item.CohortAbundance);
//                                    }
//                                }
//                            }
//                            else
//                            {
//                                foreach (int f in functionalGroups)
//                                {
//                                    foreach (var item in TempCohorts[f])
//                                    {
//                                        if (item.IndividualBodyMass <= modelInitialisation.PlanktonDispersalThreshold)
//                                            returnValue += ((item.IndividualBodyMass + item.IndividualReproductivePotentialMass) * item.CohortAbundance);
//                                    }
//                                }
//                            }
//                            break;
//
//                        case "abundance":
//                            if (traitValue != "Zooplankton (all)")
//                            {
//                                foreach (int f in functionalGroups)
//                                {
//                                    foreach (var item in TempCohorts[f])
//                                    {
//                                        returnValue += item.CohortAbundance;
//                                    }
//                                }
//                            }
//                            else
//                            {
//                                foreach (int f in functionalGroups)
//                                {
//                                    foreach (var item in TempCohorts[f])
//                                    {
//                                        if (item.IndividualBodyMass <= modelInitialisation.PlanktonDispersalThreshold)
//                                            returnValue += item.CohortAbundance;
//                                    }
//                                }
//                            }
//                            break;
//
//                        default:
//                            Debug.Fail("For cohorts, state variable name must be either 'biomass' or 'abundance'");
//                            break;
//                    }
//                    break;
//
//                case "stock":
//                    GridCellStockHandler TempStocks = InternalGrid[latCellIndex, lonCellIndex].GridCellStocks;
//
//                    switch (variableName.ToLower())
//                    {
//                        case "biomass":
//                            foreach (int f in functionalGroups)
//                            {
//                                foreach (var item in TempStocks[f])
//                                {
//                                    returnValue += item.TotalBiomass;
//                                }
//                            }
//                            break;
//                        default:
//                            Debug.Fail("For stocks, state variable name must be 'biomass'");
//                            break;
//                    }
//                    break;
//
//                default:
//                    Debug.Fail("State variable type must be either 'cohort' or 'stock'");
//                    break;
//
//            }
//
//
//            return returnValue / (InternalGrid[latCellIndex, lonCellIndex].CellEnvironment["Cell Area"][0]);
//        }
//
/** \brief
Get the mean density of a state variable for specific cells

@param variableName The name of the variable 
@param functionalGroups A vector of functional group indices to consider 
@param cellIndices List of indices of active cells in the model grid 
@param stateVariableType A string indicating the type of state variable; 'cohort' or 'stock' 
\returns Mean density of variable over whole grid
*/
//        public double StateVariableGridMeanDensity(string variableName, string traitValue, int[] functionalGroups, List<unsigned[]> cellIndices, string stateVariableType, MadingleyModelInitialisation initialisation)
//        {
//
//            double tempVal = 0;
//
//            double[,] TempStateVariable = this.GetStateVariableGridDensityPerSqKm(variableName, traitValue, functionalGroups, cellIndices, stateVariableType, initialisation);
//
//            // Loop through and sum values across a grid, excluding missing values
//            for (int ii = 0; ii < cellIndices.Count; ii++)
//            {
//                tempVal += TempStateVariable[cellIndices[ii][0], cellIndices[ii][1]];
//            }
//
//            return tempVal / cellIndices.Count;
//        }
//
//
/** \brief
Return an array of values for a single state variable over specific cells

@param variableName Variable name 
@param functionalGroups A vector of functional group indices to consider 
@param cellIndices List of indices of active cells in the model grid 
@param stateVariableType A string indicating the type of state variable; 'cohort' or 'stock' 
@return Array of state variable values for each grid cell
*/
//        public double[,] GetStateVariableGrid(string variableName, string traitValue, int[] functionalGroups, List<unsigned[]> cellIndices, string stateVariableType, MadingleyModelInitialisation initialisation)
//        {
//            double[,] TempStateVariable = new double[this.NumLatCells, this.NumLonCells];
//
//            switch (variableName.ToLower())
//            {
//                case "biomass":
//                    for (int ii = 0; ii < cellIndices.Count; ii++)
//                    {
//                        // Check whether the state variable concerns cohorts or stocks
//                        if (stateVariableType.ToLower() == "cohort")
//                        {
//                            if (traitValue != "Zooplankton")
//                            {
//                                // Check to make sure that the cell has at least one cohort
//                                if (InternalGrid[cellIndices[ii][0], cellIndices[ii][1]].GridCellCohorts != null)
//                                {
//                                    for (int nn = 0; nn < functionalGroups.Length; nn++)
//                                    {
//                                        if (InternalGrid[cellIndices[ii][0], cellIndices[ii][1]].GridCellCohorts[functionalGroups[nn]] != null)
//                                        {
//                                            foreach (Cohort item in InternalGrid[cellIndices[ii][0], cellIndices[ii][1]].GridCellCohorts[functionalGroups[nn]].ToArray())
//                                            {
//                                                    TempStateVariable[cellIndices[ii][0], cellIndices[ii][1]] += ((item.IndividualBodyMass + item.IndividualReproductivePotentialMass) * item.CohortAbundance);
//                                            }
//                                        }
//                                    }
//                                }
//                            }
//                            else
//                            {
//                                // Check to make sure that the cell has at least one cohort
//                                if (InternalGrid[cellIndices[ii][0], cellIndices[ii][1]].GridCellCohorts != null)
//                                {
//                                    for (int nn = 0; nn < functionalGroups.Length; nn++)
//                                    {
//                                        if (InternalGrid[cellIndices[ii][0], cellIndices[ii][1]].GridCellCohorts[functionalGroups[nn]] != null)
//                                        {
//                                            foreach (Cohort item in InternalGrid[cellIndices[ii][0], cellIndices[ii][1]].GridCellCohorts[functionalGroups[nn]].ToArray())
//                                            {
//                                                if (item.IndividualBodyMass <= initialisation.PlanktonDispersalThreshold)
//                                                    TempStateVariable[cellIndices[ii][0], cellIndices[ii][1]] += ((item.IndividualBodyMass + item.IndividualReproductivePotentialMass) * item.CohortAbundance);
//                                            }
//                                        }
//                                    }
//                                }
//                            }
//                        }
//                        else if (stateVariableType.ToLower() == "stock")
//                        {
//                            // Check to make sure that the cell has at least one stock
//                            if (InternalGrid[cellIndices[ii][0], cellIndices[ii][1]].GridCellStocks != null)
//                            {
//                                for (int nn = 0; nn < functionalGroups.Length; nn++)
//                                {
//                                    if (InternalGrid[cellIndices[ii][0], cellIndices[ii][1]].GridCellStocks[functionalGroups[nn]] != null)
//                                    {
//                                        foreach (Stock item in InternalGrid[cellIndices[ii][0], cellIndices[ii][1]].GridCellStocks[functionalGroups[nn]].ToArray())
//                                        {
//                                            TempStateVariable[cellIndices[ii][0], cellIndices[ii][1]] += (item.TotalBiomass);
//
//                                        }
//                                    }
//
//                                }
//                            }
//                        }
//                        else
//                        {
//                            Debug.Fail("Variable 'state variable type' must be either 'stock' 'or 'cohort'");
//                        }
//                        
//                    }
//                    break;
//                case "abundance":
//                    for (int ii = 0; ii < cellIndices.Count; ii++)
//                    {
//                        // Check whether the state variable concerns cohorts or stocks
//                        if (stateVariableType.ToLower() == "cohort")
//                        {
//                            if (traitValue != "Zooplankton")
//                            {
//                                // Check to make sure that the cell has at least one cohort
//                                if (InternalGrid[cellIndices[ii][0], cellIndices[ii][1]].GridCellCohorts != null)
//                                {
//                                    for (int nn = 0; nn < functionalGroups.Length; nn++)
//                                    {
//                                        if (InternalGrid[cellIndices[ii][0], cellIndices[ii][1]].GridCellCohorts[functionalGroups[nn]] != null)
//                                        {
//                                            foreach (Cohort item in InternalGrid[cellIndices[ii][0], cellIndices[ii][1]].GridCellCohorts[functionalGroups[nn]].ToArray())
//                                            {
//                                                TempStateVariable[cellIndices[ii][0], cellIndices[ii][1]] += item.CohortAbundance;
//                                            }
//                                        }
//                                    }
//                                }
//                            }
//                            else
//                            {
//                                // Check to make sure that the cell has at least one cohort
//                                if (InternalGrid[cellIndices[ii][0], cellIndices[ii][1]].GridCellCohorts != null)
//                                {
//                                    for (int nn = 0; nn < functionalGroups.Length; nn++)
//                                    {
//                                        if (InternalGrid[cellIndices[ii][0], cellIndices[ii][1]].GridCellCohorts[functionalGroups[nn]] != null)
//                                        {
//                                            foreach (Cohort item in InternalGrid[cellIndices[ii][0], cellIndices[ii][1]].GridCellCohorts[functionalGroups[nn]].ToArray())
//                                            {
//                                                if (item.IndividualBodyMass <= initialisation.PlanktonDispersalThreshold)
//                                                    TempStateVariable[cellIndices[ii][0], cellIndices[ii][1]] += item.CohortAbundance;
//                                            }
//                                        }
//                                    }
//                                }
//                            }
//                        }
//                        else
//                        {
//                            Debug.Fail("Currently abundance cannot be calculated for grid cell stocks");
//                        }
//                    }
//                    break;
//                default:
//                    Debug.Fail("Invalid search string passed for cohort property");
//                    break;
//            }
//
//            return TempStateVariable;
//
//        }
//        
/** \brief
Return an array of values for a single state variable over specific cells, given in densities per km^2

@param variableName Variable name 
@param functionalGroups A vector of functional group indices to consider 
@param cellIndices List of indices of active cells in the model grid 
@param stateVariableType A string indicating the type of state variable; 'cohort' or 'stock' 
@return Array of state variable values for each grid cell
*/
//        public double[,] GetStateVariableGridDensityPerSqKm(string variableName, string traitValue, int[] functionalGroups, List<unsigned[]> cellIndices, string stateVariableType, MadingleyModelInitialisation initialisation)
//        {
//            double[,] TempStateVariable = new double[this.NumLatCells, this.NumLonCells];
//            double CellArea;
//
//            TempStateVariable = this.GetStateVariableGrid(variableName, traitValue, functionalGroups, cellIndices, stateVariableType, initialisation);
//
//            for (int ii = 0; ii < cellIndices.Count; ii++)
//            {
//                CellArea = GetCellEnvironment(cellIndices[ii][0], cellIndices[ii][1])["Cell Area"][0];
//                TempStateVariable[cellIndices[ii][0], cellIndices[ii][1]] /= CellArea;
//            }
//
//            return TempStateVariable;
//        }
//
//
/** \brief
Return an array of log(values + 1) for a state variable for particular functional groups over specific cells. State variable (currently only biomass or abundance) must be >= 0 in all grid cells

@param variableName The name of the variable 
@param functionalGroups A vector of functional group indices to consider 
@param cellIndices List of indices of active cells in the model grid 
@param stateVariableType A string indicating the type of state variable; 'cohort' or 'stock' 
@return Array of log(state variable values +1 ) for each grid cell
*/
//        public double[,] GetStateVariableGridLog(string variableName, string traitValue, int[] functionalGroups, List<unsigned[]> cellIndices, string stateVariableType, MadingleyModelInitialisation initialisation)
//        {
//
//            double[,] TempStateVariable = new double[this.NumLatCells, this.NumLonCells];
//
//            TempStateVariable = this.GetStateVariableGrid(variableName, traitValue, functionalGroups, cellIndices, stateVariableType, initialisation);
//            
//            for (int ii = 0; ii < cellIndices.Count; ii++)
//            {
//                TempStateVariable[cellIndices[ii][0], cellIndices[ii][1]] = Math.Log(TempStateVariable[cellIndices[ii][0], cellIndices[ii][1]]+1);
//            }
//
//            return TempStateVariable;
//        }
//
//
/** \brief
Return an array of log(values + 1) for a state variable for particular functional groups over specific cells. State variable (currently only biomass or abundance) must be >= 0 in all grid cells

@param variableName The name of the variable 
@param functionalGroups A vector of functional group indices to consider 
@param cellIndices List of indices of active cells in the model grid 
@param stateVariableType A string indicating the type of state variable; 'cohort' or 'stock' 
@return Array of log(state variable values +1 ) for each grid cell
*/
//        public double[,] GetStateVariableGridLogDensityPerSqKm(string variableName, string traitValue, int[] functionalGroups, List<unsigned[]> cellIndices, string stateVariableType, MadingleyModelInitialisation initialisation)
//        {
//
//            double[,] TempStateVariable = new double[this.NumLatCells, this.NumLonCells];
//            double CellArea;
//
//            TempStateVariable = this.GetStateVariableGrid(variableName, traitValue, functionalGroups, cellIndices, stateVariableType, initialisation);
//
//            for (int ii = 0; ii < cellIndices.Count; ii++)
//            {
//                CellArea = GetCellEnvironment(cellIndices[ii][0], cellIndices[ii][1])["Cell Area"][0];
//                TempStateVariable[cellIndices[ii][0], cellIndices[ii][1]] /= CellArea;
//                TempStateVariable[cellIndices[ii][0], cellIndices[ii][1]] = Math.Log(TempStateVariable[cellIndices[ii][0], cellIndices[ii][1]]+1);
//            }
//
//            return TempStateVariable;
//
//        }
//
/** \brief
Returns, for a given longitude, the appropriate longitude index in the grid
ASSUMES THAT LONGITUDES IN THE MODEL GRID OBJECT REFER TO LOWER LEFT CORNERS!!!

@param myLon Longitude, in degrees 
@return longitude index in the model grid
*/
//        public unsigned GetLonIndex(double myLon)
//        {
//            Debug.Assert((myLon >= _MinLongitude && myLon < _MaxLongitude), "Error: latitude out of range");
//
//            return (unsigned)Math.Floor((myLon - _MinLongitude) / _LonCellSize);
//
//        }
//
/** \brief
Return the longitude of a cell at a particular lon. index

@param cellLonIndex The longitudinal index (i.e. row) of the cell 
@return Returns the longitude of the bottom of the cell, in degrees
*/
//        public double GetCellLongitude(unsigned cellLonIndex)
//        {
//            Debug.Assert((cellLonIndex <= (_NumLonCells - 1)), "Error: Cell index out of range when trying to find the longitude for a particular cell");
//
//            double TempLongitude = double.MaxValue;
//
//            for (int ii = 0; ii < _NumLatCells; ii++)
//            {
//                if (InternalGrid[ii, cellLonIndex] != null)
//                    TempLongitude = InternalGrid[ii, cellLonIndex].Longitude;
//            }
//
//            Debug.Assert(TempLongitude != double.MaxValue, "Error trying to find cell longitude - no grid cells have been initialised for this latitude index: " + cellLonIndex.ToString());
//
//            return TempLongitude;
//        }
//
/** \brief
Return the latitude of a cell at a particular lat. index

@param cellLatIndex The latitudinal index (i.e. row) of the cell 
@return Returns the latitude of the bottom of the cell, in degrees
*/
//        public double GetCellLatitude(unsigned cellLatIndex)
//        {
//            Debug.Assert((cellLatIndex <= (_NumLatCells - 1)), "Error: Cell index out of range when trying to find the latitude for a particular cell");
//
//            double TempLatitude = double.MaxValue;
//
//            for (int jj = 0; jj < _NumLonCells ; jj++)
//            {
//                if (InternalGrid[cellLatIndex, jj] != null)
//                {
//                    TempLatitude = InternalGrid[cellLatIndex, jj].Latitude;
//                    break;
//                }
//            }
//
//            Debug.Assert(TempLatitude != double.MaxValue, "Error trying to find cell latitude - no grid cells have been initialised for this latitude index: " + cellLatIndex.ToString());
//
//            return TempLatitude;
//
//        }
//
//
/** \brief
Returns, for a given latitude, the appropriate latitude index in the grid
ASSUMES THAT LATITUDES IN THE MODEL GRID OBJECT REFER TO LOWER LEFT CORNERS!!!

@param myLat Latitude, in degrees 
@return latitude index in the model grid
*/
//        public unsigned GetLatIndex(double myLat)
//        {
//
//            Debug.Assert((myLat >= _MinLatitude && myLat < _MaxLatitude), "Error: latitude out of range");
//
//            return (unsigned)Math.Floor((myLat - _MinLatitude) / _LatCellSize);
//
//        }
//
/** \brief
A method to return the values for all environmental data layers for a particular grid cell
@param cellLatIndex Latitude index of grid cell 
@param cellLonIndex Longitude index of grid cell 
@return A sorted list containing environmental data layer names and values
*/
map<string, vector<double> > GetCellEnvironment(unsigned cellLatIndex, unsigned cellLonIndex)
        {
//            return InternalGrid[cellLatIndex, cellLonIndex].CellEnvironment;
        }
//
/** \brief
A method to return delta values for the specified delta type in a particular grid cell
@param deltaType The delta type to return 
@param cellLatIndex Latitude index of grid cell 
@param cellLonIndex Longitude index of grid cell 
@return A sorted list containing deltas
*/
//        public Dictionary<string, double> GetCellDeltas(string deltaType, unsigned cellLatIndex, unsigned cellLonIndex)
//        {
//            return InternalGrid[cellLatIndex, cellLonIndex].Deltas[deltaType];
//        }
//
/** \brief
A method to return all delta values in a particular grid cell
@param cellLatIndex Latitude index of grid cell 
@param cellLonIndex Longitude index of grid cell 
@return A sorted list of sorted lists containing deltas
*/
//        public Dictionary<string, Dictionary<string, double>> GetCellDeltas(unsigned cellLatIndex, unsigned cellLonIndex)
//        {
//            return InternalGrid[cellLatIndex, cellLonIndex].Deltas;
//        }
//
//        
/** \brief
Get a grid of values for an environmental data layer
@param enviroVariable  The name of the environmental data layer 
@param timeInterval The desired time interval within the environmental variable (i.e. 0 if it is a yearly variable
or the month index - 0=Jan, 1=Feb etc. - for monthly variables) 
@return The values in each grid cell
*/
//        public double[,] GetEnviroGrid(string enviroVariable,unsigned timeInterval)
//        {
//            // Check to see if environmental variable exists
//            for (int ii = 0; ii < _NumLatCells; ii++)
//            {
//                for (int jj = 0; jj < _NumLonCells; jj++)
//                {
//                    if(InternalGrid[ii,jj] != null)
//                        Debug.Assert(InternalGrid[ii, jj].CellEnvironment.ContainsKey(enviroVariable), "Environmental variable not found when running GetEnviroGrid");
//                }
//            }
//
//            double[,] outputData = new double[_NumLatCells, _NumLonCells];
//
//            for (int ii = 0; ii < _NumLatCells; ii+=GridCellRarefaction)
//            {
//                for (int jj = 0; jj < _NumLonCells; jj+=GridCellRarefaction)
//                {
//                    outputData[ii, jj] = InternalGrid[ii, jj].CellEnvironment[enviroVariable][timeInterval];
//                }
//            }
//
//            return outputData;
//        }
//
/** \brief
Get a grid of values for an environmental data layer in specific cells
@param enviroVariable The name of the environmental data layer to return 
@param timeInterval The desired time interval for which to get data (i.e. 0 if it is a yearly variable
or the month index - 0=Jan, 1=Feb etc. - for monthly variables) 
@param cellIndices List of active cells in the model grid 
@return The values in each grid cell
*/
//        public double[,] GetEnviroGrid(string enviroVariable, unsigned timeInterval, List<unsigned[]> cellIndices)
//        {
//            // Check to see if environmental variable exists
//            for (int ii = 0; ii < cellIndices.Count; ii++)
//            {
//                if (InternalGrid[cellIndices[ii][0], cellIndices[ii][1]] != null)
//                        Debug.Assert(InternalGrid[cellIndices[ii][0], cellIndices[ii][1]].CellEnvironment.ContainsKey(enviroVariable), 
//                            "Environmental variable not found when running GetEnviroGrid");
//                
//            }
//
//            // Create grid to hold the data to return
//            double[,] outputData = new double[_NumLatCells, _NumLonCells];
//
//            for (int ii = 0; ii < cellIndices.Count; ii++)
//            {
//                outputData[cellIndices[ii][0], cellIndices[ii][1]] = InternalGrid[cellIndices[ii][0], cellIndices[ii][1]].CellEnvironment
//                    [enviroVariable][timeInterval];
//            }
//
//            return outputData;
//        }
//
/** \brief
Return the total over the whole grid for an environmental variable
@param enviroVariable The environmental variable 
@param timeInterval The desired time interval within the environmental variable (i.e. 0 if it is a yearly variable
or the month index - 0=Jan, 1=Feb etc. - for monthly variables) 
@return The total of the variable over the whole grid
*/
//        public double GetEnviroGridTotal(string enviroVariable, unsigned timeInterval)
//        {
//            double[,] enviroGrid = GetEnviroGrid(enviroVariable,timeInterval);
//            double enviroTotal = 0.0;
//
//            for (int ii = 0; ii < _NumLatCells; ii+=GridCellRarefaction)
//            {
//                for (int jj = 0; jj < _NumLonCells; jj+=GridCellRarefaction)
//                {
//                    enviroTotal += enviroGrid[ii, jj];
//                }
//            }
//
//            return enviroTotal;
//        }
//
/** \brief
Return the sum of an environmental variable over specific cells
@param enviroVariable The environmental variable 
@param timeInterval The desired time interval within the environmental variable (i.e. 0 if it is a yearly variable
or the month index - 0=Jan, 1=Feb etc. - for monthly variables) 
@param cellIndices List of active cells in the model grid 
@return The total of the variable over the whole grid
*/
//        public double GetEnviroGridTotal(string enviroVariable, unsigned timeInterval, List<unsigned[]> cellIndices)
//        {
//            double[,] enviroGrid = GetEnviroGrid(enviroVariable,timeInterval, cellIndices);
//            double enviroTotal = 0.0;
//
//            for (int ii = 0; ii < cellIndices.Count; ii++)
//            {
//                enviroTotal += enviroGrid[cellIndices[ii][0], cellIndices[ii][1]];
//            }
//
//            return enviroTotal;
//        }
//
/** \brief
Check to see if the top perimeter of the cell is traversable for dispersal (i.e. is from the same realm)
@param latCell The latitudinal cell index 
@param lonCell The longitudinal cell index 
@param gridCellRealm The grid cell realm 
*/
//        private void CheckTopPerimeterTraversable(unsigned latCell, unsigned lonCell, double gridCellRealm)
//        {
//            // Check to see if top perimeter is traversable
//            if (InternalGrid[latCell + 1, lonCell].CellEnvironment["Realm"][0] == gridCellRealm)
//            {
//                // Add the cell above to the list of cells that are dispersable to
//                CellsForDispersal[latCell, lonCell].Add(new unsigned[2] { (latCell + 1), (lonCell) });
//
//                // Also add it to the directional list
//                CellsForDispersalDirection[latCell, lonCell].Add(1);
//            }
//        }
//
/** \brief
Check to see if the top right perimeter of the cell is traversable for dispersal (i.e. is from the same realm)
@param latCell The latitudinal cell index 
@param lonCell The longitudinal cell index 
@param lonCellToGoTo The index of the cell to go to (needs to take into account grid wrapping) 
@param gridCellRealm The grid cell realm 
*/
//        private void CheckTopRightPerimeterTraversable(unsigned latCell, unsigned lonCell, unsigned lonCellToGoTo, double gridCellRealm)
//        {
//
//            // Check to see if right perimeter is traversable
//            if (InternalGrid[latCell + 1, lonCellToGoTo].CellEnvironment["Realm"][0] == gridCellRealm)
//            {
//                // Add the cell above to the list of cells that are dispersable to
//                CellsForDispersal[latCell, lonCell].Add(new unsigned[2] { (latCell + 1), (lonCellToGoTo) });
//                
//                // Also add it to the directional list
//                CellsForDispersalDirection[latCell, lonCell].Add(2);
//            }
//        }
//
/** \brief
Check to see if the right perimeter of the cell is traversable for dispersal (i.e. is from the same realm)
@param latCell The latitudinal cell index 
@param lonCell The longitudinal cell index 
@param lonCellToGoTo The index of the cell to go to (needs to take into account grid wrapping) 
@param gridCellRealm The grid cell realm 
*/
//        private void CheckRightPerimeterTraversable(unsigned latCell, unsigned lonCell, unsigned lonCellToGoTo, double gridCellRealm)
//        {
//            // Check to see if right perimeter is traversable
//            if (InternalGrid[latCell, lonCellToGoTo].CellEnvironment["Realm"][0] == gridCellRealm)
//            {
//                // Add the cell above to the list of cells that are dispersable to
//                CellsForDispersal[latCell, lonCell].Add(new unsigned[2] { (latCell), (lonCellToGoTo) });
//
//                // Also add it to the directional list
//                CellsForDispersalDirection[latCell, lonCell].Add(3);
//
//            }
//        }
//
//
/** \brief
Check to see if the bottom right perimeter of the cell is traversable for dispersal (i.e. is from the same realm)
@param latCell The latitudinal cell index 
@param lonCell The longitudinal cell index 
@param lonCellToGoTo The index of the cell to go to (needs to take into account grid wrapping) 
@param gridCellRealm The grid cell realm 
*/
//        private void CheckBottomRightPerimeterTraversable(unsigned latCell, unsigned lonCell, unsigned lonCellToGoTo, double gridCellRealm)
//        {
//            // Check to see if bottom right perimeter is traversable
//            if (InternalGrid[latCell - 1, lonCellToGoTo].CellEnvironment["Realm"][0] == gridCellRealm)
//            {
//                // Add the cell above to the list of cells that are dispersable to
//                CellsForDispersal[latCell, lonCell].Add(new unsigned[2] { (latCell - 1), (lonCellToGoTo) });
//                
//                // Also add it to the directional list
//                CellsForDispersalDirection[latCell, lonCell].Add(4);
//            }
//        }
//
//
/** \brief
Check to see if the right perimeter of the cell is traversable for dispersal (i.e. is from the same realm)
@param latCell The latitudinal cell index 
@param lonCell The longitudinal cell index 
@param gridCellRealm The grid cell realm 
*/
//        private void CheckBottomPerimeterTraversable(unsigned latCell, unsigned lonCell, double gridCellRealm)
//        {
//            // Check to see if top perimeter is traversable
//            if (InternalGrid[latCell - 1, lonCell].CellEnvironment["Realm"][0] == gridCellRealm)
//            {
//                // Add the cell above to the list of cells that are dispersable to
//                CellsForDispersal[latCell, lonCell].Add(new unsigned[2] { (latCell - 1), (lonCell) });
//
//                // Also add it to the directional list
//                CellsForDispersalDirection[latCell, lonCell].Add(5);
//            
//            }
//        }
//
//
/** \brief
Check to see if the bottom left perimeter of the cell is traversable for dispersal (i.e. is from the same realm)
@param latCell The latitudinal cell index 
@param lonCell The longitudinal cell index 
@param lonCellToGoTo The index of the cell to go to (needs to take into account grid wrapping) 
@param gridCellRealm The grid cell realm 
*/
//        private void CheckBottomLeftPerimeterTraversable(unsigned latCell, unsigned lonCell, unsigned lonCellToGoTo, double gridCellRealm)
//        {
//            // Check to see if bottom right perimeter is traversable
//            if (InternalGrid[latCell - 1, lonCellToGoTo].CellEnvironment["Realm"][0] == gridCellRealm)
//            {
//                // Add the cell above to the list of cells that are dispersable to
//                CellsForDispersal[latCell, lonCell].Add(new unsigned[2] { (latCell - 1), (lonCellToGoTo) });
//                
//                // Also add it to the directional list
//                CellsForDispersalDirection[latCell, lonCell].Add(6);
//
//            }
//        }
//
//
/** \brief
Check to see if the left perimeter of the cell is traversable for dispersal (i.e. is from the same realm)
@param latCell The latitudinal cell index 
@param lonCell The longitudinal cell index 
@param lonCellToGoTo The index of the cell to go to (needs to take into account grid wrapping) 
@param gridCellRealm The grid cell realm 
*/
//        private void CheckLeftPerimeterTraversable(unsigned latCell, unsigned lonCell, unsigned lonCellToGoTo, double gridCellRealm)
//        {
//            // Check to see if left perimeter is traversable
//            if (InternalGrid[latCell, lonCellToGoTo].CellEnvironment["Realm"][0] == gridCellRealm)
//            {
//                // Add the cell above to the list of cells that are dispersable to
//                CellsForDispersal[latCell, lonCell].Add(new unsigned[2] { (latCell), (lonCellToGoTo) });
//
//                // Also add it to the directional list
//                CellsForDispersalDirection[latCell, lonCell].Add(7);
//
//            }
//        }
//
//
/** \brief
Check to see if the top left perimeter of the cell is traversable for dispersal (i.e. is from the same realm)
@param latCell The latitudinal cell index 
@param lonCell The longitudinal cell index 
@param lonCellToGoTo The index of the cell to go to (needs to take into account grid wrapping) 
@param gridCellRealm The grid cell realm 
*/
//        private void CheckTopLeftPerimeterTraversable(unsigned latCell, unsigned lonCell, unsigned lonCellToGoTo, double gridCellRealm)
//        {
//            // Check to see if bottom right perimeter is traversable
//            if (InternalGrid[latCell + 1, lonCellToGoTo].CellEnvironment["Realm"][0] == gridCellRealm)
//            {
//                // Add the cell above to the list of cells that are dispersable to
//                CellsForDispersal[latCell, lonCell].Add(new unsigned[2] { (latCell + 1), (lonCellToGoTo) });
//
//                // Also add it to the directional list
//                CellsForDispersalDirection[latCell, lonCell].Add(8);
//
//            }
//        }
// 
//        // Currently assumes that the grid does not run from -90 to 90 (in which case there would be transfer at top and bottom latitude)
//        // Also needs checking to see if it works with a sub-grid
/** \brief
Calculate the dispersable perimeter lengths of each of the grid cells
*/
//        private void CalculatePerimeterLengthsAndCellsDispersableTo()
//        {
//            int counter = 1;
//            // Loop through grid cells
//            for (unsigned ii = 0; ii < _NumLatCells; ii++)
//			{
//                // Bottom of the grid
//			    if (ii == 0)
//                {
//                    // Loop through the longitude indices of each cell
//                    for (unsigned jj = 0; jj < _NumLonCells; jj++)
//			        {
//                        // Get the realm of the cell (i.e. whether it is land or sea)
//                        double GridCellRealm = InternalGrid[ii,jj].CellEnvironment["Realm"][0];
//                        if ((GridCellRealm != 1.0) && (GridCellRealm != 2.0))
//                        {
//                            Console.Write("\r{0} cells classified as neither land nor sea",counter);
//                            counter++;
//                            break;
//                        }
//
//                        // Check to see if we are at the left-most edge
//			            if (jj == 0)
//                        {
//                            // Are we on a grid that spans the globe?
//                            if ((_MaxLongitude - _MinLongitude) > 359.9)
//                            {
//                                // Check to see if the top perimeter is dispersable
//                                CheckTopPerimeterTraversable(ii, jj, GridCellRealm);
//
//                                // Check to see if the top right perimeter is dispersable
//                                CheckTopRightPerimeterTraversable(ii, jj, jj + 1, GridCellRealm);
//                                
//                                // Check to see if the right perimeter is dispersable
//                                CheckRightPerimeterTraversable(ii, jj, jj + 1, GridCellRealm);
//
//                                // Check to see if the left perimeter is dispersable
//                                CheckLeftPerimeterTraversable(ii, jj, _NumLonCells - 1, GridCellRealm);
//
//                                // Check to see if the top left perimeter is dispersable
//                                CheckTopLeftPerimeterTraversable(ii, jj, _NumLonCells - 1, GridCellRealm);                  
//                            }
//
//                            // Otherwise, we are simply on a non-wrappable boundary. 
//                            // Assumes that we have a closed system on this boundary and that organisms cannot disperse through it
//                            else
//                            {
//                                // Check to see if the top perimeter is traversable
//                                CheckTopPerimeterTraversable(ii, jj, GridCellRealm);
//                                
//                                // Check to see if the top right perimeter is dispersable
//                                CheckTopRightPerimeterTraversable(ii, jj, jj + 1, GridCellRealm);
//                                
//                                // Check to see if the right perimeter is dispersable
//                                CheckRightPerimeterTraversable(ii, jj, jj + 1, GridCellRealm);
//                            }
//                        }
//                        // Check to see if we are at the right-most edge
//                        else if (jj == (_NumLonCells - 1))
//                        {
//                            // Are we on a grid that spans the globe?
//                            if ((_MaxLongitude - _MinLongitude) > 359.9)
//                            {
//                                // Check to see if the top perimeter is traversable
//                                CheckTopPerimeterTraversable(ii, jj, GridCellRealm);
//
//                                // Check to see if the top right perimeter is dispersable
//                                CheckTopRightPerimeterTraversable(ii, jj, 0, GridCellRealm);
//                                
//                                // Check to see if the right perimeter is dispersable
//                                CheckRightPerimeterTraversable(ii, jj, 0, GridCellRealm);
//
//                                // Check to see if the left perimeter is dispersable
//                                CheckLeftPerimeterTraversable(ii, jj, jj - 1, GridCellRealm);
//
//                                // Check to see if the top left perimeter is dispersable
//                                CheckTopLeftPerimeterTraversable(ii, jj, jj - 1, GridCellRealm);
//                            }
//                            // Otherwise, we are simply on a non-wrappable boundary. 
//                            // Assumes that we have a closed system on this boundary and that organisms cannot disperse through it
//                            else
//                            {
//                                // Check to see if the top perimeter is traversable
//                                CheckTopPerimeterTraversable(ii, jj, GridCellRealm);
//
//                                // Check to see if the left perimeter is dispersable
//                                CheckLeftPerimeterTraversable(ii, jj, jj - 1, GridCellRealm);
//
//                                // Check to see if the top left perimeter is dispersable
//                                CheckTopLeftPerimeterTraversable(ii, jj, jj - 1, GridCellRealm);
//                            }
//                        }
//
//                        // Otherwise internal in the grid longitudinally
//                        else
//                        {
//                            // Check to see if the top perimeter is traversable
//                            CheckTopPerimeterTraversable(ii, jj, GridCellRealm);
//
//                            // Check to see if the top right perimeter is dispersable
//                            CheckTopRightPerimeterTraversable(ii, jj, jj + 1, GridCellRealm);
//                                
//                            // Check to see if the right perimeter is dispersable
//                            CheckRightPerimeterTraversable(ii, jj, jj + 1, GridCellRealm);
//
//                            // Check to see if the left perimeter is dispersable
//                            CheckLeftPerimeterTraversable(ii, jj, jj - 1, GridCellRealm);
//
//                            // Check to see if the top left perimeter is dispersable
//                            CheckTopLeftPerimeterTraversable(ii, jj, jj - 1, GridCellRealm);
//                        }
//			        }
//                }
//
//                // Top of the grid
//                else if (ii == (_NumLatCells -1))
//                {
//                    // Loop through the longitude indices of each cell
//                    for (unsigned jj = 0; jj < _NumLonCells; jj++)
//			        {
//                        // Get the realm of the cell (i.e. whether it is land or sea)
//                        double GridCellRealm = InternalGrid[ii,jj].CellEnvironment["Realm"][0];
//                        if ((GridCellRealm != 1.0) && (GridCellRealm != 2.0))
//                        {
//                            Console.Write("\r{0} cells classified as neither land nor sea", counter);
//                            counter++;
//                            break;
//                        }
//                    
//                        // Check to see if we are at the left-most edge
//                        if (jj == 0)
//                        {
//                              // Are we on a grid that spans the globe?
//                                if ((_MaxLongitude - _MinLongitude) > 359.9)
//                                {
//                                    // Check to see if the right perimeter is dispersable
//                                    CheckRightPerimeterTraversable(ii, jj, jj + 1, GridCellRealm);
//
//                                    // Check to see if the bottom right perimeter is dispersable
//                                    CheckBottomRightPerimeterTraversable(ii, jj, jj + 1, GridCellRealm);
//                                
//                                    // Check to see if the bottom perimeter is dispersable
//                                    CheckBottomPerimeterTraversable(ii, jj, GridCellRealm);
//
//                                    // Check to see if the bottom left perimeter is dispersable
//                                    CheckBottomLeftPerimeterTraversable(ii, jj, _NumLonCells - 1, GridCellRealm);
//
//                                    // Check to see if the left perimeter is dispersable
//                                    CheckLeftPerimeterTraversable(ii, jj, _NumLonCells - 1, GridCellRealm);
//                                }
//                                // Otherwise, we are simply on a non-wrappable boundary. 
//                                // Assumes that we have a closed system on this boundary and that organisms cannot disperse through it
//                                else
//                                {
//                                    // Check to see if the right perimeter is dispersable
//                                    CheckRightPerimeterTraversable(ii, jj, jj + 1, GridCellRealm);
//                                    
//                                    // Check to see if the bottom right perimeter is dispersable
//                                    CheckBottomRightPerimeterTraversable(ii, jj, jj + 1, GridCellRealm);
//
//                                    // Check to see if the bottom perimeter is dispersable
//                                    CheckBottomPerimeterTraversable(ii, jj, GridCellRealm);
//                                }
//                        }
//                        // Check to see if we are at the right-most edge
//                        else if (jj == (_NumLonCells - 1))
//                        {
//                            // Are we on a grid that spans the globe?
//                            if ((_MaxLongitude - _MinLongitude) > 359.9)
//                            {
//                                // Check to see if the right perimeter is dispersable
//                                CheckRightPerimeterTraversable(ii, jj, 0, GridCellRealm);
//                                
//                                // Check to see if the bottom right perimeter is dispersable
//                                CheckBottomRightPerimeterTraversable(ii, jj, 0, GridCellRealm);
//
//                                // Check to see if the bottom perimeter is dispersable
//                                CheckBottomPerimeterTraversable(ii, jj, GridCellRealm);
//                                
//                                // Check to see if the bottom left perimeter is dispersable
//                                CheckBottomLeftPerimeterTraversable(ii, jj, jj - 1, GridCellRealm);
//
//                                // Check to see if the left perimeter is dispersable
//                                CheckLeftPerimeterTraversable(ii, jj, jj - 1, GridCellRealm);
//                            }
//
//                            // Otherwise, we are simply on a non-wrappable boundary. 
//                            // Assumes that we have a closed system on this boundary and that organisms cannot disperse through it
//                            else
//                            {
//                                // Check to see if the bottom perimeter is dispersable
//                                CheckBottomPerimeterTraversable(ii, jj, GridCellRealm);
//
//                                // Check to see if the bottom left perimeter is dispersable
//                                CheckBottomLeftPerimeterTraversable(ii, jj, jj - 1, GridCellRealm);
//
//                                // Check to see if the left perimeter is dispersable
//                                CheckLeftPerimeterTraversable(ii, jj, jj - 1, GridCellRealm);
//                            }
//
//                        }
//                        // Otherwise, internal in the grid longitudinally
//                        else
//                        {
//                            // Check to see if the right perimeter is dispersable
//                            CheckRightPerimeterTraversable(ii, jj, jj + 1, GridCellRealm);
//
//                            // Check to see if the bottom right perimeter is dispersable
//                            CheckBottomRightPerimeterTraversable(ii, jj, jj + 1, GridCellRealm);
//
//                            // Check to see if the bottom perimeter is dispersable
//                            CheckBottomPerimeterTraversable(ii, jj, GridCellRealm);
//
//                            // Check to see if the bottom left perimeter is dispersable
//                            CheckBottomLeftPerimeterTraversable(ii, jj, jj - 1, GridCellRealm);
//
//                            // Check to see if the left perimeter is dispersable
//                            CheckLeftPerimeterTraversable(ii, jj, jj - 1, GridCellRealm);
//                        }
//                    }
//
//                }
//                // Otherwise internal latitudinally
//                else
//                {
//                    // Loop through the longitude indices of each cell
//                    for (unsigned jj = 0; jj < _NumLonCells; jj++)
//                    {
//                        // Get the realm of the cell (i.e. whether it is land or sea)
//                        double GridCellRealm = InternalGrid[ii, jj].CellEnvironment["Realm"][0];
//                        if ((GridCellRealm != 1.0) && (GridCellRealm != 2.0))
//                        {
//                            Console.Write("\r{0} cells classified as neither land nor sea", counter);
//                            counter++;
//                            break;
//                        }
//
//                        // Check to see if we are at the left-most edge
//                        if (jj == 0)
//                        {
//                            // Are we on a grid that spans the globe?
//                            if ((_MaxLongitude - _MinLongitude) > 359.9)
//                            {
//                                // Check to see if the top perimeter is dispersable
//                                CheckTopPerimeterTraversable(ii, jj, GridCellRealm);
//
//                                // Check to see if the top right perimeter is dispersable
//                                CheckTopRightPerimeterTraversable(ii, jj, jj + 1, GridCellRealm);
//
//                                // Check to see if the right perimeter is dispersable
//                                CheckRightPerimeterTraversable(ii, jj, jj + 1, GridCellRealm);
//
//                                // Check to see if the bottom right perimeter is dispersable
//                                CheckBottomRightPerimeterTraversable(ii, jj, jj + 1, GridCellRealm);
//
//                                // Check to see if the bottom perimeter is dispersable
//                                CheckBottomPerimeterTraversable(ii, jj, GridCellRealm);
//                                
//                                // Check to see if the bottom left perimeter is dispersable
//                                CheckBottomLeftPerimeterTraversable(ii, jj, _NumLonCells - 1, GridCellRealm);
//
//                                // Check to see if the left perimeter is dispersable
//                                CheckLeftPerimeterTraversable(ii, jj, _NumLonCells - 1, GridCellRealm);
//
//                                // Check to see if the top left perimeter is dispersable
//                                CheckTopLeftPerimeterTraversable(ii, jj, _NumLonCells - 1, GridCellRealm);
//                            }
//                            // Otherwise, we are simply on a non-wrappable boundary. 
//                            // Assumes that we have a closed system on this boundary and that organisms cannot disperse through it
//                            else
//                            {
//                                // Check to see if the top perimeter is dispersable
//                                CheckTopPerimeterTraversable(ii, jj, GridCellRealm);
//
//                                // Check to see if the top right perimeter is dispersable
//                                CheckTopRightPerimeterTraversable(ii, jj, jj + 1, GridCellRealm);
//
//                                // Check to see if the right perimeter is dispersable
//                                CheckRightPerimeterTraversable(ii, jj, jj + 1, GridCellRealm);
//                                
//                                // Check to see if the bottom right perimeter is dispersable
//                                CheckBottomRightPerimeterTraversable(ii, jj, jj + 1, GridCellRealm);
//
//                                // Check to see if the bottom perimeter is dispersable
//                                CheckBottomPerimeterTraversable(ii, jj, GridCellRealm);
//                            }
//                        }
//                        // Check to see if we are at the rightmost edge
//                        else if (jj == (_NumLonCells - 1))
//                        {
//                            // Are we on a grid that spans the globe?
//                            if ((_MaxLongitude - _MinLongitude) > 359.9)
//                            {
//                                // Check to see if the top perimeter is dispersable
//                                CheckTopPerimeterTraversable(ii, jj, GridCellRealm);
//                                
//                                // Check to see if the top right perimeter is dispersable
//                                CheckTopRightPerimeterTraversable(ii, jj, 0, GridCellRealm);
//
//                                // Check to see if the right perimeter is dispersable
//                                CheckRightPerimeterTraversable(ii, jj, 0, GridCellRealm);
//                                
//                                // Check to see if the bottom right perimeter is dispersable
//                                CheckBottomRightPerimeterTraversable(ii, jj, 0, GridCellRealm);
//
//                                // Check to see if the bottom perimeter is dispersable
//                                CheckBottomPerimeterTraversable(ii, jj, GridCellRealm);
//                                
//                                // Check to see if the bottom left perimeter is dispersable
//                                CheckBottomLeftPerimeterTraversable(ii, jj, jj - 1, GridCellRealm);
//
//                                // Check to see if the left perimeter is dispersable
//                                CheckLeftPerimeterTraversable(ii, jj, jj - 1, GridCellRealm);
//                                
//                                // Check to see if the top left perimeter is dispersable
//                                CheckTopLeftPerimeterTraversable(ii, jj, jj - 1, GridCellRealm);
//                            }
//                            else
//                            {
//                                // Check to see if the top perimeter is dispersable
//                                CheckTopPerimeterTraversable(ii, jj, GridCellRealm);
//
//                                // Check to see if the bottom perimeter is dispersable
//                                CheckBottomPerimeterTraversable(ii, jj, GridCellRealm);
//
//                                // Check to see if the bottom left perimeter is dispersable
//                                CheckBottomLeftPerimeterTraversable(ii, jj, jj - 1, GridCellRealm);
//
//                                // Check to see if the left perimeter is dispersable
//                                CheckLeftPerimeterTraversable(ii, jj, jj - 1, GridCellRealm);
//
//                                // Check to see if the top left perimeter is dispersable
//                                CheckTopLeftPerimeterTraversable(ii, jj, jj - 1, GridCellRealm);
//                            }
//                        }
//                        // Otherwise internal in the grid both latitudinally and longitudinally - the easiest case
//                        else
//                        {
//                            // Check to see if the top perimeter is dispersable
//                            CheckTopPerimeterTraversable(ii, jj, GridCellRealm);
//
//                            // Check to see if the top right perimeter is dispersable
//                            CheckTopRightPerimeterTraversable(ii, jj, jj + 1, GridCellRealm);
//
//                            // Check to see if the right perimeter is dispersable
//                            CheckRightPerimeterTraversable(ii, jj, jj + 1, GridCellRealm);
//
//                            // Check to see if the bottom right perimeter is dispersable
//                            CheckBottomRightPerimeterTraversable(ii, jj, jj + 1, GridCellRealm);
//
//                            // Check to see if the bottom perimeter is dispersable
//                            CheckBottomPerimeterTraversable(ii, jj, GridCellRealm);
//
//                            // Check to see if the bottom left perimeter is dispersable
//                            CheckBottomLeftPerimeterTraversable(ii, jj, jj - 1, GridCellRealm);
//
//                            // Check to see if the left perimeter is dispersable
//                            CheckLeftPerimeterTraversable(ii, jj, jj - 1, GridCellRealm);
//
//                            // Check to see if the top left perimeter is dispersable
//                            CheckTopLeftPerimeterTraversable(ii, jj, jj - 1, GridCellRealm);
//                        }
//                    }
//                 }
//			}
//            Console.WriteLine("\n");
//
//        }
//
/** \brief
Given a grid cell from where a cohort is dispersing, select at random a grid cell for it to disperse to from those that exist within the 
same realm

@param fromCellLatIndex The latitudinal index of the cell from which the cohort is dispersing 
@param fromCellLonIndex The longitudinal index of the cell from which the cohort is dispersing 
@return what?
*/
vector<unsigned> GetRandomGridCellToDisperseTo(unsigned fromCellLatIndex, unsigned fromCellLonIndex)
        {
//            // Select a cell at random
//            int CellPickedAtRandom = (int)Math.Floor(RandomNumberGenerator.GetUniform() * CellsForDispersal[fromCellLatIndex, fromCellLonIndex].Count);
//
//            // Return the coordinates of that cell
//            return CellsForDispersal[fromCellLatIndex, fromCellLonIndex][CellPickedAtRandom];
        }
//
/** \brief
Get the longitudinal and latitudinal indices of the cell that lies to the north of the focal grid cell, if a viable cell to disperse to

@param fromCellLatIndex The latitudinal index of the focal grid cell 
@param fromCellLonIndex The longitudinal index of the focal grid cell 
@return The longitudinal and latitudinal cell indcies of the cell that lies to the north of the focal grid cell
*/
vector<unsigned> CheckDispersalNorth(unsigned fromCellLatIndex, unsigned fromCellLonIndex)
        {
//            unsigned[] NorthCell = new unsigned[2] { 9999999, 9999999 };
//
//            for (int ii = 0; ii < CellsForDispersalDirection[fromCellLatIndex, fromCellLonIndex].Count; ii++)
//            {
//                if (CellsForDispersalDirection[fromCellLatIndex, fromCellLonIndex][ii] == 1)
//                {
//                    NorthCell = CellsForDispersal[fromCellLatIndex, fromCellLonIndex][ii];
//                    break;
//                }
//            }
//
//            return NorthCell;
        }
//
/** \brief
Get the longitudinal and latitudinal indices of the cell that lies to the east of the focal grid cell, if a viable cell to disperse to
@param fromCellLatIndex The latitudinal index of the focal grid cell 
@param fromCellLonIndex The longitudinal index of the focal grid cell 
@return The longitudinal and latitudinal cell indcies of the cell that lies to the east of the focal grid cell

*/
vector<unsigned> CheckDispersalEast(unsigned fromCellLatIndex, unsigned fromCellLonIndex)
        {
//            unsigned[] EastCell = new unsigned[2] {9999999,9999999};
//
//            for (int ii = 0; ii < CellsForDispersalDirection[fromCellLatIndex, fromCellLonIndex].Count; ii++)
//            {
//                if (CellsForDispersalDirection[fromCellLatIndex, fromCellLonIndex][ii] == 3)
//                {
//                    EastCell = CellsForDispersal[fromCellLatIndex, fromCellLonIndex][ii];
//                    break;
//                }
//            }
//
//            return EastCell;
        }
//
/** \brief
Get the longitudinal and latitudinal indices of the cell that lies to the south of the focal grid cell, if a viable cell to disperse to
@param fromCellLatIndex The latitudinal index of the focal grid cell 
@param fromCellLonIndex The longitudinal index of the focal grid cell 
@return The longitudinal and latitudinal cell indcies of the cell that lies to the south of the focal grid cell
*/
vector<unsigned>  CheckDispersalSouth(unsigned fromCellLatIndex, unsigned fromCellLonIndex)
        {
//            unsigned[] SouthCell = new unsigned[2] { 9999999, 9999999 };
//
//            for (int ii = 0; ii < CellsForDispersalDirection[fromCellLatIndex, fromCellLonIndex].Count; ii++)
//            {
//                if (CellsForDispersalDirection[fromCellLatIndex, fromCellLonIndex][ii] == 5)
//                {
//                    SouthCell = CellsForDispersal[fromCellLatIndex, fromCellLonIndex][ii];
//                }
//            }
//
//            return SouthCell;
        }
//
/** \brief
Get the longitudinal and latitudinal indices of the cell that lies to the west of the focal grid cell, if a viable cell to disperse to
@param fromCellLatIndex The latitudinal index of the focal grid cell 
@param fromCellLonIndex The longitudinal index of the focal grid cell 
@return The longitudinal and latitudinal cell indcies of the cell that lies to the west of the focal grid cell
*/
vector<unsigned> CheckDispersalWest(unsigned fromCellLatIndex, unsigned fromCellLonIndex)
        {
//            unsigned[] WestCell = new unsigned[2] { 9999999, 9999999 };
//
//            for (int ii = 0; ii < CellsForDispersalDirection[fromCellLatIndex, fromCellLonIndex].Count; ii++)
//            {
//                if (CellsForDispersalDirection[fromCellLatIndex, fromCellLonIndex][ii] == 7)
//                {
//                    WestCell = CellsForDispersal[fromCellLatIndex, fromCellLonIndex][ii];
//                }
//            }
//
//            return WestCell;
        }
//
/** \brief
Get the longitudinal and latitudinal indices of the cell that lies to the northeast of the focal grid cell, if a viable cell to disperse to

@param fromCellLatIndex The latitudinal index of the focal grid cell 
@param fromCellLonIndex The longitudinal index of the focal grid cell 
@return The longitudinal and latitudinal cell indcies of the cell that lies to the northeast of the focal grid cell
*/
vector<unsigned>  CheckDispersalNorthEast(unsigned fromCellLatIndex, unsigned fromCellLonIndex)
        {
//            unsigned[] NECell = new unsigned[2] { 9999999, 9999999 };
//
//            for (int ii = 0; ii < CellsForDispersalDirection[fromCellLatIndex, fromCellLonIndex].Count; ii++)
//            {
//                if (CellsForDispersalDirection[fromCellLatIndex, fromCellLonIndex][ii] == 2)
//                {
//                    NECell = CellsForDispersal[fromCellLatIndex, fromCellLonIndex][ii];
//                }
//            }
//
//            return NECell;
        }
//
/** \brief
Get the longitudinal and latitudinal indices of the cell that lies to the southeast of the focal grid cell, if a viable cell to disperse to

@param fromCellLatIndex The latitudinal index of the focal grid cell 
@param fromCellLonIndex The longitudinal index of the focal grid cell 
@return The longitudinal and latitudinal cell indcies of the cell that lies to the southeast of the focal grid cell
*/
vector<unsigned> CheckDispersalSouthEast(unsigned fromCellLatIndex, unsigned fromCellLonIndex)
        {
//            unsigned[] SECell = new unsigned[2] { 9999999, 9999999 };
//
//            for (int ii = 0; ii < CellsForDispersalDirection[fromCellLatIndex, fromCellLonIndex].Count; ii++)
//            {
//                if (CellsForDispersalDirection[fromCellLatIndex, fromCellLonIndex][ii] == 4)
//                {
//                    SECell = CellsForDispersal[fromCellLatIndex, fromCellLonIndex][ii];
//                }
//            }
//
//            return SECell;
        }
//
/** \brief
Get the longitudinal and latitudinal indices of the cell that lies to the southwest of the focal grid cell, if a viable cell to disperse to

@param fromCellLatIndex The latitudinal index of the focal grid cell 
@param fromCellLonIndex The longitudinal index of the focal grid cell 
@return The longitudinal and latitudinal cell indcies of the cell that lies to the southwest of the focal grid cell
*/
vector<unsigned> CheckDispersalSouthWest(unsigned fromCellLatIndex, unsigned fromCellLonIndex)
        {
//            unsigned[] SWCell = new unsigned[2] { 9999999, 9999999 };
//
//            for (int ii = 0; ii < CellsForDispersalDirection[fromCellLatIndex, fromCellLonIndex].Count; ii++)
//            {
//                if (CellsForDispersalDirection[fromCellLatIndex, fromCellLonIndex][ii] == 6)
//                {
//                    SWCell = CellsForDispersal[fromCellLatIndex, fromCellLonIndex][ii];
//                }
//            }
//
//            return SWCell;
        }
//
/** \brief
Get the longitudinal and latitudinal indices of the cell that lies to the northwest of focal grid cell, if a viable cell to disperse to

@param fromCellLatIndex The latitudinal index of the focal grid cell 
@param fromCellLonIndex The longitudinal index of the focal grid cell 
@return The longitudinal and latitudinal cell indcies of the cell that lies to the northwest of the focal grid cell
*/
vector<unsigned> CheckDispersalNorthWest(unsigned fromCellLatIndex, unsigned fromCellLonIndex)
        {
//            unsigned[] NWCell = new unsigned[2] { 9999999, 9999999 };
//
//            for (int ii = 0; ii < CellsForDispersalDirection[fromCellLatIndex, fromCellLonIndex].Count; ii++)
//            {
//                if (CellsForDispersalDirection[fromCellLatIndex, fromCellLonIndex][ii] == 8)
//                {
//                    NWCell = CellsForDispersal[fromCellLatIndex, fromCellLonIndex][ii];
//                }
//            }
//
//            return NWCell;
//        }
//
    }
};
#endif