#ifndef ENVIRODATA_H
#define ENVIRODATA_H

//    /// <summary>
//    /// Imports environmental data from ASCII and NetCDF files
//    /// </summary>
//    /// <todoT>No error-trapping as yet</todoT>
//    /// <todoT>Rewrite to use the ArraySDSConvert class</todoT>
//    /// <todoD>Need  to go through code and rewrite e.g. change method to overloaded to prevent passing variable name and file name for ESRI grids</todoD>
//    /// <remarks>Currently assumes that cells are evenly spaced in latitude and longitude</remarks>
class EnviroData
    {
    public:
            //----------------------------------------------------------------------------------------------
    //Variables
    //----------------------------------------------------------------------------------------------
/** \brief
Number of latitudinal cells
*/
         unsigned NumLats;
//
/** \brief
Number of longitudinal cells
*/
          unsigned NumLons;

/** brief Number of time intervals encompassed by the environmental variable */
        unsigned NumTimes;
/** \brief
Latitude of the bottom edge of the sothernmost grid cell
*/
//        private double _LatMin;
/** \brief
Get latitude of the bottom edge of the sothernmost grid cell
*/
//        public double LatMin { get { return _LatMin; } }
//
/** \brief
Latitude of the left edge of the most western grid cell
*/
//        private double _LonMin;
/** \brief
Get latitude of the left edge of the most western grid cell
*/
//        public double LonMin { get { return _LonMin; } }
//
/** \brief
Value used to denote missing data for this environmental variable
*/
//        private double _MissingValue;
/** \brief
Get value used to denote missing data for this environmental variable
*/
//        public double MissingValue { get { return _MissingValue; } }
//
/** \brief
Latitudinal distance between adjacent cells
*/
//        private double _LatStep;
/** \brief
Get latitudinal distance between adjacent cells
*/
//        public double LatStep { get { return _LatStep; } }
//
/** \brief
Longitudinal distance between adjacent cells
*/
//        private double _LonStep;
/** \brief
Get longitudinal distance between adjacent cells
*/
//        public double LonStep { get { return _LonStep; } }
//
/** \brief
List of arrays of values of the environmental variable
*/
vector<vector<vector<double>>> DataArray;//time by lon by lat
/** \brief
Get list of arrays of values of the environmental variable
*/
//        public List<double[,]> DataArray { get { return _DataArray; } }
//
/** \brief
Vector of latitudes of the bottom edges of grid cells
*/
//        private double[] _Lats;
/** \brief
Get vector of latitudes of the bottom edges of grid cells
*/
//        public double[] Lats { get { return _Lats; } }
//
/** \brief
Vector of longitudes of the left edges of grid cells
*/
//        private double[] _Lons;
/** \brief
Get vector of longitudes of the left edges of grid cells
*/
//        public double[] Lons { get { return _Lons; } }
//
/** \brief
Vector containing values of the time dimension of the environmental variable
*/
//        private double[] _Times;
/** \brief
Get vector containing values of the time dimension of the environmental variable
*/
//        public double[] Times { get { return _Times; } }
//
/** \brief
The string required to read the file with the environmental data
*/
//        private string _ReadFileString;
/** \brief
Get the string required to read the file with the environmental data
*/
//        public string ReadFileString { get { return _ReadFileString; } }
//
/** \brief
The units of the environmental variable
*/
//        private string _Units;
/** \brief
Gets the units of the environmental variable
*/
//        public string Units
//        { get { return _Units; } }
//
/** \brief
Tracks the number of environmental data layers opened
*/
//        private uint _NumEnviroLayers;
/** \brief
Returns the number of environmental data layers opened
*/
//        public uint NumEnviroLayers { get { return _NumEnviroLayers; } }
//
/** \brief
Instance of the class to perform general functions
*/
//        private UtilityFunctions Utilities;
//
    //----------------------------------------------------------------------------------------------
    //Methods
    //----------------------------------------------------------------------------------------------
EnviroData(){;}
/** \brief
Constructor for EnviroData

@param fileName Filename (including extension) 
@param dataName The name of the variable that contains the data within the specified file 
@param dataType Type of data, nc = NetCDF, ascii = ESRI ASCII) 
@param dataResolution The temporal resolution of the environmental variable 
@param units The units of the data 
<todo>Check whether lat/lon or 0/1 are fixed for all NetCDFs</todo>
<todo>CHECK IF DIMENSIONS HAVE TO BE THE SAME FOR ALL VARIABLES IN A NETCDF AND HOW TO EXTRACT DIMENSIONS FOR A SINGLE VARIABLE IF NECESSARY</todo>
<todo>Write code to check for equal cell sizes in NetCDFs</todo>
*/
EnviroData(string fileName, string dataName, string dataType, string dataResolution, string units)
       {
    //dummy data for initial testing
    int NumLonCells=1;int NumLatCells=1;double value=0;
    if(dataName=="SST")value=293;
    if(dataName=="AWC")value=0.1;
    if(dataName=="NPP")value=0.1;
    if(dataName=="land_sea_mask")value=0;
    if (dataName=="frost")value=0;
    vector<vector< double> > TempStateVariable(NumLatCells);for (auto &t : TempStateVariable )t.resize(NumLonCells);
    for (unsigned i=0;i<NumLatCells;i++)for (unsigned j=0;j<NumLonCells;j++)TempStateVariable[i][j]=value;
    for (int i=0;i<12;i++)DataArray.push_back(TempStateVariable);
    NumTimes=12;
        }

/** \brief
A method to extract the area weighted value of an environmental variable from the envirodata cells overlapped by the cell specified by lat and lon

@param lat Bottom latitude of cell to get value from 
@param lon Leftmost longitude of cell to get value from 
@param timeInterval The time interval to get the value from (i.e. the month, or 0 for yearly variables) 
@param missingValue Boolean to indicate whether the returned value is a missing value 
@param latCellSize The latitudinal size of cells in the model grid 
@param lonCellSize The longitudinal size of cells in the model grid 
@return The area weighted value of an environmental variable from the envirodata cells overlapped by the cell specified by lat and lon
*/
double GetValue(double lat, double lon, unsigned timeInterval, bool missingValue, double latCellSize, double lonCellSize)
       {
        return DataArray[0][0][0];
        missingValue=false;
       }

};
#endif