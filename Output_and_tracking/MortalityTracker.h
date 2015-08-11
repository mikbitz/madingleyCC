#ifndef MORTALITYTRACKER_H
#define MORTALITYTRACKER_H
/** \file MortalityTracker.h
 * \brief the MortalityTracker header file
 */






//

//
//namespace Madingley
//{
/** \brief
//    /// Track results associated with mortality
//    /// </summary>
//     class MortalityTracker
//    {
/** \brief
Name of the file to write data on mortality to
*/
//        string MortalityFilename;
//
/** \brief
A streamwriter instance to output data on mortality
*/
//        private StreamWriter MortalityWriter;
//
/** \brief
Synchronized version of the streamwriter to output mortality data
*/
//        private TextWriter SyncMortalityWriter;
//
/** \brief
List of mortality events for each cohort (keyed by cohort ID)
*/
//        private map<string,List<string>> MortalityList;
//
//
/** \brief
An instance of the simple random number generator
*/
//        private NonStaticRNG RandomNumberGenerator = new NonStaticRNG();
//
/** \brief
Set up properties of the mortality tracker
*/
@param numTimeSteps The total number of time steps for this simulation 
@param numLats The number of latitudinal cells in the model grid 
@param numLons The number of longitudinal cells in the model grid 
@param cellIndices List of indices of active cells in the model grid 
@param mortalityFilename The filename to write data on mortality to 
@param outputFileSuffix The suffix to apply to all output files from this model run 
@param outputPath The path to write all output files to 
//         MortalityTracker(unsigned numTimeSteps,
//            unsigned numLats, unsigned numLons,
//            List<unsigned[]> cellIndices,
//            string mortalityFilename,
//            string outputFileSuffix,
//            string outputPath, int cellIndex)
//        {
//            MortalityFilename = mortalityFilename;
//
//            // Initialise streamwriter to output mortality of cohorts
//            MortalityWriter = new StreamWriter(outputPath + MortalityFilename + outputFileSuffix + "_Cell" + cellIndex + ".txt");
//            // Create a threadsafe textwriter to write outputs to the NewCohortWriter stream
//            SyncMortalityWriter = TextWriter.Synchronized(MortalityWriter);
//            SyncMortalityWriter.WriteLine("Latitude\tLongitude\tbirth_step\ttime_step\tcurrent mass\tadult mass\tfunctional group\tcohort id\tnumber died\tmortality source");
//
//            MortalityList = new map<string,List<string>>();
//
//        
//        }
//
/** \brief
Record a mortality event associated with a cohort to memory
*/
@param latIndex The latitudinal index of the current grid cell 
@param lonIndex The longitudinal index of the current grid cell 
@param birthStep The time step that the cohort came into existence 
@param timestep The current time step 
@param currentMass The current body mass of individuals in the cohort with dying individuals 
@param adultMass The adult body mass of individuals in the cohort with dying individuals 
@param functionalGroup The index of the functional group that the cohort belongs to 
@param cohortID The unique ID of the cohort 
@param numberDied The number of individuals  
@param mortalitySource  
//         void RecordMortality(unsigned latIndex, unsigned lonIndex,
//            unsigned birthStep, unsigned timestep, double currentMass,double adultMass, unsigned functionalGroup, 
//            unsigned cohortID, double numberDied, string mortalitySource)
//        {
//            // Write the time step and the abundance of the new cohort to the output file for diagnostic purposes
//            string newline = Convert.ToString(latIndex + "\t" + lonIndex + "\t" + birthStep + "\t" + timestep + "\t" + currentMass + "\t" + adultMass + "\t" +
//            functionalGroup + "\t" + cohortID + "\t" + numberDied + "\t" + mortalitySource);
//
//            string CohortIDString = Convert.ToString(cohortID);
//
//            if (MortalityList.ContainsKey(CohortIDString))
//            {
//                MortalityList[CohortIDString].Add(newline);
//            }
//            else
//            {
//                MortalityList.Add(CohortIDString, new List<string>() { newline });
//            }
//
//            /*
//            string StringCID = Convert.ToString(cohortID);
//
//            lock (MortalityList.SyncRoot)
//            {
//                if (MortalityList.ContainsKey(StringCID))
//                {
//                    //Copy the current contents for this cohortID to a variable
//                    List<object> objects = ((IEnumerable)MortalityList[StringCID]).Cast<object>().ToList();
//                    List<string> strings = (from o in objects select o.ToString()).ToList();
//
//                    //Add to the list the newline
//                    strings.Add(newline);
//                    //Remove the current MortalityList entry for cohortID
//                    MortalityList.Remove(StringCID);
//                    //Add the new lines to MortalityList for this cohortID
//                    MortalityList.Add(StringCID, strings);
//
//                }
//                else
//                {
//                    if (RandomNumberGenerator.GetUniform() > 0.99)
//                    {
//                        //Add a new MortalityList entry for cohortID
//                        List<string> strings = new List<string>();
//                        strings.Add(newline);
//
//                        MortalityList.Add(StringCID, strings);
//                    }
//                }
//            }
//            */
//        }
//
/** \brief
Output the mortality profile of a cohort becoming extinct
*/
@param cohortID The ID of the cohort becoming extinct 
//         void OutputMortalityProfile(unsigned cohortID)
//        {
//            
//            string CohortIDString = Convert.ToString(cohortID);
//
//            if (MortalityList.ContainsKey(CohortIDString))
//            {
//                if (RandomNumberGenerator.GetUniform() > 0.95)
//                {
//                    var Lines = ((IEnumerable)MortalityList[CohortIDString]).Cast<object>().ToList();
//                    foreach (var Line in Lines)
//                    {
//                        SyncMortalityWriter.WriteLine(Line);
//                    }
//                }
//                MortalityList.Remove(CohortIDString);
//            }
//
//            /*
//            if (MortalityList.ContainsKey(StringCID))
//            {
//                lock (MortalityList.SyncRoot)
//                {
//                    var Lines = ((IEnumerable)MortalityList[StringCID]).Cast<object>().ToList();
//                    // Loop over mortality events for the specified cohort and write these to the output file
//                    foreach (var Line in Lines)
//                    {
//                        SyncMortalityWriter.WriteLine(Line);
//                    }
//
//                    // Remove the cohort from the list of mortality profiles
//                    MortalityList.Remove(Convert.ToString(cohortID));
//                }
//            }
//            */
//        }
//
//
//    }
//}
#endif
