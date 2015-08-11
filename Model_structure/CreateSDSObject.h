#ifndef CREATESDSOBJECT_H
#define CREATESDSOBJECT_H
/** \file CreateSDSObject.h
 * \brief the CreateSDSObject header file
 */






//




//
//namespace Madingley
//{
/** \brief
//    /// Class for creating Scientific Dataset objects
//    /// </summary>
//     class CreateSDSObject
//    {
/** \brief
Create an SDS object in memory
*/
@param shared Boolean indicating whether a shared dataset is required 
@return The new dataset object*/
//         DataSet CreateSDSInMemory(bool shared)
//        {
//            // Create and SDS object
//            DataSet internalSDS;
//
//            // If a shared dataset has been specified, then open a shared dataset object, otherwise open an ordinary dataset object
//            if (shared)
//            {
//                internalSDS = SharedDataSet.Open("msds:memory2");
//            }
//            else
//            {
//                internalSDS = DataSet.Open("msds:memory2");
//            }
//
//            // Disable auto commit
//            internalSDS.IsAutocommitEnabled = false;
//
//            // Return the SDS object
//            return internalSDS;
//        }
//
/** \brief
Create an SDS object as an output file
*/
@param sdsType The type of output file to create, currently must be NetCDF 
@param sdsName The name to assign to the output file 
@param outputPath The path to the output folder 
@return The new dataset object*/
//         DataSet CreateSDS(string sdsType, string sdsName, string outputPath)
//        {
//            // Check that the user has not specified an SDS object of type memory
//            if (sdsType == "Memory")
//                Debug.Fail("Error: you do not need to specify a file name for SDS objects of type 'memory'");
//            
//            // Check that the output file does not already exist
//            if (sdsType == "netCDF")
//            {
//                string filePath = outputPath + sdsName + ".nc";
//                if (System.IO.File.Exists(filePath))
//                    Debug.Fail("Error: SDS object already exists");
//            }
//
//            // If the output type had been selected as NetCDF, then create the object, otherwise throw an error for now
//            if (sdsType == "netCDF")
//            {
//                // Create the URI for the SDS object to be created
//                string tempString = "msds:nc?file="+outputPath + sdsName + ".nc&openMode=create";
//                // Create an SDS object
//                DataSet internalSDS = DataSet.Open(tempString);
//                // Disable auto commit
//                internalSDS.IsAutocommitEnabled = false;
//                // Return the new SDS object
//                return internalSDS;
//            }
//            else
//            {
//                // Throw an error
//                Debug.Fail("Error: specified type not supported. Only 'netCDF' is supported at present");
//                // Nonsense return
//                DataSet internalSDS = DataSet.Open("nonsense");
//                return internalSDS;
//            }
//        }
//    }
//}
#endif
