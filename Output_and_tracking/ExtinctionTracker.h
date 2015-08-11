#ifndef EXTINCTIONTRACKER_H
#define EXTINCTIONTRACKER_H
/** \file ExtinctionTracker.h
 * \brief the ExtinctionTracker header file
 */





//

//
//namespace Madingley
//{
/** \brief
//    /// Tracks results associated with cohort extinction
//    /// </summary>
//     class ExtinctionTracker
//    {
//        string ExtinctionFilename;
//
//        private StreamWriter ExtinctionWriter;
//
//        private TextWriter SyncedExtinctionWriter;
//
/** \brief
Constructor for the eating tracker: sets up output file
*/
@param extinctionFilename The filename for the output file 
@param outputPath The path to the output directory 
@param outputFilesSuffix The suffix to be applied to all outputs from this model simulation 
@param cellIndex The index of the current cell within the list of cells in this simulation 
//         ExtinctionTracker(string extinctionFilename, string outputPath, string outputFilesSuffix, int cellIndex)
//        {
//            ExtinctionFilename = extinctionFilename;
//
//            // Initialise streamwriter to output properties and ids of extinct cohorts
//            ExtinctionWriter = new StreamWriter(outputPath + extinctionFilename + outputFilesSuffix + "_Cell" + cellIndex + ".txt");
//            // Create a threadsafe textwriter to write outputs to the ExtinctionWriter stream
//            SyncedExtinctionWriter = TextWriter.Synchronized(ExtinctionWriter);
//            SyncedExtinctionWriter.WriteLine("Latitude\tLongitude\ttime_step\tmerged\tcohortID");
//
//        }
//
/** \brief
Record the extinction of a cohort in the output file
*/
@param latIndex The latitudinal index of the current grid cell 
@param lonIndex The longitudinal index of the current grid cell 
@param currentTimeStep The current model time step 
@param merged Whether the cohort going extinct has ever been merged with another cohort 
@param cohortID The ID of the cohort going extinct 
//         void RecordExtinction(unsigned latIndex, unsigned lonIndex,unsigned currentTimeStep,bool merged,List<unsigned> cohortID)
//        {
//            string newline = Convert.ToString(latIndex) + '\t' + Convert.ToString(lonIndex) + '\t' +
//                Convert.ToString(currentTimeStep) + '\t' + Convert.ToString(merged) + '\t' +
//                Convert.ToString(cohortID[0]);
//
//            SyncedExtinctionWriter.WriteLine(newline);
//        }
//
//
//
//    }
//}
#endif
