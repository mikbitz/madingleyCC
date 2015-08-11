#ifndef NPPTRACKER_H
#define NPPTRACKER_H
/** \file NPPTracker.h
 * \brief the NPPTracker header file
 */





//

//
//namespace Madingley
//{
/** \brief
//    /// Tracks primary productivity
//    /// </summary>
//     class NPPTracker
//    {
/** \brief
The filename for the NPP data
*/
//        string NPPFileName;
//
/** \brief
Stream-writer to output NPP data
*/
//        private StreamWriter NPPWriter;
//
/** \brief
Thread-safe text-writer to output NPP data
*/
//        private TextWriter SyncedNPPWriter;
//
/** \brief
Set up the tracker for outputting NPP data to file
*/
@param nppFilename The name of the file to write information on NPP to 
@param outputPath The file path to write all outputs to 
@param outputFilesSuffix The suffix to apply to output files from this simulation 
//         NPPTracker(string nppFilename, string outputPath, string outputFilesSuffix)
//        {
//            NPPFileName = nppFilename;
//
//            // Initialise stream-writers to output NPP data
//            NPPWriter = new StreamWriter(outputPath + NPPFileName + outputFilesSuffix + ".txt");
//            SyncedNPPWriter = TextWriter.Synchronized(NPPWriter);
//            SyncedNPPWriter.WriteLine("Latitude\tLongitude\ttime_step\tcell_area\ttotal_cell_productivity_g_per_month");
//
//        }
//
/** \brief
Record the total primary productivity in the current cell in the current time step
*/
@param latIndex The latitudinal index of the current cell 
@param lonIndex The longitudinal index of the current cell 
@param timeStep The current model time step 
@param cellArea The area of the current grid cell 
@param cellNPP The total primary productivity in the cell this time step 
//         void RecordNPP(unsigned latIndex, unsigned lonIndex, unsigned timeStep, double cellArea, double cellNPP)
//        {
//            SyncedNPPWriter.WriteLine(Convert.ToString(latIndex) + '\t' + Convert.ToString(lonIndex) + '\t' + Convert.ToString(timeStep) +
//                '\t' + Convert.ToString(cellArea) + '\t' + Convert.ToString(cellNPP));
//
//        }
//
/** \brief
Close the streams for writing NPP data
*/
//         void CloseStreams()
//        {
//            SyncedNPPWriter.Dispose();
//            NPPWriter.Dispose();
//        }
//        
//    }
//}
#endif
