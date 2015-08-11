#ifndef METABOLISMTRACKER_H
#define METABOLISMTRACKER_H
/** \file MetabolismTracker.h
 * \brief the MetabolismTracker header file
 */






//
//namespace Madingley
//{
/** \brief
//    /// Tracks variables associated with cohort metabolism
//    /// </summary>
//     class MetabolismTracker
//    {
/** \brief
File to write metabolism data to
*/
//        string MetabolismFilename;
//
/** \brief
A streamwriter to write out data on metabolism
*/
//        private StreamWriter MetabolismWriter;
//        private TextWriter SyncMetabolismWriter;
//
/** \brief
Set up the metabolism tracker
*/
@param metabolismFilename Name of the metabolism tracker file to write to 
@param outputPath The path to the folder in which the metabolism tracker data file will be stored 
@param outputFilesSuffix A suffix for the filename in the case that there is more than one scenario 
@param cellIndex The index of the current cell within the list of all cells in this simulation 
//         MetabolismTracker(string metabolismFilename, string outputPath, string outputFilesSuffix, int cellIndex)
//        {
//            MetabolismFilename = metabolismFilename;
//
//            MetabolismWriter = new StreamWriter(outputPath + MetabolismFilename + outputFilesSuffix + "_Cell" + cellIndex + ".txt");
//            SyncMetabolismWriter = TextWriter.Synchronized(MetabolismWriter);
//            SyncMetabolismWriter.WriteLine("Latitude\tLongitude\ttime_step\tCurrent_body_mass\tfunctional_group\tAmbient_temp\tMetabolic_mass_loss");
//        }
//
/** \brief
Record the metabolic loss of individuals in a cohort
*/
@param latIndex The latitudinal index of the current grid cell 
@param lonIndex The longitudinal index of the current grid cell 
@param timeStep The current time step 
@param currentBodyMass The current body mass of individuals in the cohort 
@param functionalGroup The index of the functional group that the cohort belongs to 
@param temperature The ambient temperature this cohort is experiencing 
@param metabolicLoss The metabolic loss of this cohort in this time step 
//         void RecordMetabolism(unsigned latIndex, unsigned lonIndex, unsigned timeStep, double currentBodyMass, int functionalGroup, double temperature, double metabolicLoss)
//        {
//            SyncMetabolismWriter.WriteLine(Convert.ToString(latIndex) + "\t" +
//                                            Convert.ToString(lonIndex) + "\t" +
//                                            Convert.ToString(timeStep) + "\t" +
//                                            Convert.ToString(currentBodyMass) + "\t" +
//                                            Convert.ToString(functionalGroup) + "\t" +
//                                            Convert.ToString(temperature) + "\t" +
//                                            Convert.ToString(metabolicLoss));
//        }
//
/** \brief
Closes streams for writing metabolism data
*/
//         void CloseStreams()
//        {
//
//            SyncMetabolismWriter.Dispose();
//            MetabolismWriter.Dispose();
//        }
//    }
//}
#endif
