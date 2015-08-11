#ifndef GLOBALNPPTRACKER_H
#define GLOBALNPPTRACKER_H
/** \file GlobalNPPTracker.h
 * \brief the GlobalNPPTracker header file
 */







//
//namespace Madingley
//{
//     class GlobalNPPTracker
//    {
//        //An array to hold the global NPP fields
//        double[,] NPP;
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
A dataset to store the NPP outputs to file
*/
//        private DataSet NPPOutput;
//
//        private int _NumLats;
//        private int _NumLons;
//
//         GlobalNPPTracker(string outputPath,int numLats, int numLons, float[] lats, float[] lons, float latCellSize,float lonCellSize, int numTimeSteps)
//        {
//            _NumLats = numLats;
//            _NumLons = numLons;
//
//            // Initialise the data converter
//            DataConverter = new ArraySDSConvert();
//
//            // Initialise the SDS object creator
//            SDSCreator = new CreateSDSObject();
//
//            // Create an SDS object to hold total abundance and biomass data
//            NPPOutput = SDSCreator.CreateSDS("netCDF", "NPPOutput", outputPath);
//
//            // Create vector to hold the values of the time dimension
//            float[] TimeSteps = new float[numTimeSteps];
//
//
//
//            // Fill other values from 0 (this will hold outputs during the model run)
//            for (int i = 0; i < numTimeSteps; i++)
//            {
//                TimeSteps[i] = i;
//            }
//
//            // Declare vectors for geographical dimension data
//            float[] outLats = new float[numLats];
//            float[] outLons = new float[numLons];
//
//            // Populate the dimension variable vectors with cell centre latitude and longitudes
//            for (int i = 0; i < numLats; i++)
//            {
//                outLats[i] = lats[i] + (latCellSize / 2);
//            }
//
//            for (int jj = 0; jj < numLons; jj++)
//            {
//                outLons[jj] = lons[jj] + (lonCellSize / 2);
//            }
//
//
//            // Add output variables that are dimensioned geographically and temporally to grid output file
//            string[] GeographicalDimensions = { "Latitude", "Longitude", "Time step" };
//            DataConverter.AddVariable(NPPOutput, "NPP", 3, GeographicalDimensions, -9999.0, outLats, outLons, TimeSteps);
//
//
//            NPP = new double[numLats, numLons];
//            for (int ii = 0; ii < numLats; ii++)
//            {
//                for (int jj = 0; jj < numLons; jj++)
//                {
//                    NPP[ii,jj] = -9999.0;
//                }
//                
//            }
//
//        }
//
/** \brief
Add the NPP value for this grid cell
*/
@param latIndex The latitude index of the grid cell 
@param lonIndex The longitude index of the grid cell 
@param val The NPP value to be recorded 
//         void RecordNPPValue(unsigned latIndex,unsigned lonIndex, double val)
//        {
//            NPP[latIndex, lonIndex] = val;
//        }
//
/** \brief
Add the filled NPP grid the memory dataset ready to be written to file
*/
@param t  
//         void StoreNPPGrid(unsigned t)
//        {
//            DataConverter.Array2DToSDS3D(NPP, "NPP", new string[] { "Latitude", "Longitude", "Time step" },
//                                        (int)t, 0, NPPOutput);
//
//            for (int ii = 0; ii < _NumLats; ii++)
//            {
//                for (int jj = 0; jj < _NumLons; jj++)
//                {
//                    NPP[ii, jj] = -9999.0;
//                }
//
//            }
//        }
//
//         void CloseNPPFile()
//        {
//            NPPOutput.Dispose();
//        }
//
//    }
//}
#endif
