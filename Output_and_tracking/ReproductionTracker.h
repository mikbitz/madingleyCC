#ifndef REPRODUCTIONTRACKER_H
#define REPRODUCTIONTRACKER_H
/** \file ReproductionTracker.h
 * \brief the ReproductionTracker header file
 */








//
//namespace Madingley
//{
/** \brief
//    /// Tracks results associated with the reproduction process
//    /// </summary>
//     class ReproductionTracker
//    {
/** \brief
File to write data on newly produced cohorts to
*/
//        string NewCohortsFilename;
//
/** \brief
File to write data on maturity of cohorts to
*/
//        string MaturityFilename;
//
/** \brief
A streamwriter instance for outputting data on newly produced cohorts
*/
//        private StreamWriter NewCohortWriter;
/** \brief
Synchronized version of the streamwriter for outputting data on newly produced cohorts
*/
//        private TextWriter SyncNewCohortWriter;
/** \brief
A streamwriter instance for outputting data on maturity of cohorts
*/
//        private StreamWriter MaturityWriter;
/** \brief
A synchronized version of the streamwriter for outuputting data on the maturity of cohorts
*/
//        private TextWriter SyncMaturityWriter;
//
/** \brief
Sets up properties of the reproduction tracker
*/
@param numTimeSteps The total number of timesteps for this simulation 
@param numLats The number of latitudes in the model grid 
@param numLons The number of longitudes in the model grid 
@param cellIndices List of indices of active cells in the model grid 
@param newCohortsFilename The filename to write information about new cohorts to 
@param maturityFilename The filename to write information about cohorts reaching maturity 
@param outputFileSuffix The suffix to apply to all output files from this model run 
@param outputPath The path to write all output files to 
//         ReproductionTracker(unsigned numTimeSteps,
//            unsigned numLats, unsigned numLons, 
//            List<unsigned[]> cellIndices, 
//            string newCohortsFilename, 
//            string maturityFilename,
//            string outputFileSuffix,
//            string outputPath, int cellIndex)
//        {
//
//            NewCohortsFilename = newCohortsFilename;
//            MaturityFilename = maturityFilename;
//            
//
//            // Initialise streamwriter to output abundance of newly produced cohorts to a text file
//            NewCohortWriter = new StreamWriter(outputPath + newCohortsFilename + outputFileSuffix + "_Cell" + cellIndex + ".txt");
//            // Create a threadsafe textwriter to write outputs to the NewCohortWriter stream
//            SyncNewCohortWriter = TextWriter.Synchronized(NewCohortWriter);
//            SyncNewCohortWriter.WriteLine("Latitude\tLongitude\ttime_step\tabundance\tfunctional group\tadult mass\tparent cohort IDs\toffspring cohort ID");
//
//            MaturityWriter = new StreamWriter(outputPath + maturityFilename + outputFileSuffix + "_Cell" + cellIndex + ".txt");
//            // Create a threadsafe textwriter to write outputs to the Maturity stream
//            SyncMaturityWriter = TextWriter.Synchronized(MaturityWriter);
//            SyncMaturityWriter.WriteLine("Latitude\tLongitude\ttime_step\tbirth_step\tjuvenile Mass\tadult mass\tfunctional group");
//            
//        }
//
/** \brief
Records information about new cohorts spawned in the model
*/
@param latIndex The latitude index of the grid cell in which the cohort was spawned 
@param lonIndex The longitude index of the grid cell in which the cohort was spawned 
@param timestep The model timestep in which the spawning happened 
@param offspringCohortAbundance The abundance of the offspring cohort 
@param parentCohortAdultMass The adult mass of the parent cohort 
@param functionalGroup The functional group of the offspring cohort 
@param parentCohortIDs The cohort IDs associated with the parent cohort 
@param offspringCohortID The cohort ID used for the new offspring cohort 
//         void RecordNewCohort(unsigned latIndex, unsigned lonIndex, unsigned timestep, double offspringCohortAbundance, double parentCohortAdultMass, 
//            int functionalGroup, List<unsigned> parentCohortIDs,unsigned offspringCohortID)
//        {
//            double[] NewCohortRecords = new double[3];
//            NewCohortRecords[0] = offspringCohortAbundance;
//            NewCohortRecords[1] = parentCohortAdultMass;
//            NewCohortRecords[2] = (double)functionalGroup;
//
//            string AllCohortIDs = Convert.ToString(parentCohortIDs[0]);
//            if (parentCohortIDs.Count > 1)
//            {
//                for (int i = 1; i < parentCohortIDs.Count; i++)
//                {
//                    AllCohortIDs = AllCohortIDs + "; " + Convert.ToString(parentCohortIDs[i]);
//                }
//            }
//
//            // Write the time step and the abundance of the new cohort to the output file for diagnostic purposes
//            string newline = Convert.ToString(latIndex) + '\t' + Convert.ToString(lonIndex) + '\t' +
//                Convert.ToString(timestep) + '\t' + Convert.ToString(offspringCohortAbundance) + '\t' +
//                Convert.ToString(functionalGroup) + '\t' + Convert.ToString(parentCohortAdultMass) + '\t' + AllCohortIDs +
//                '\t' + Convert.ToString(offspringCohortID);
//            SyncNewCohortWriter.WriteLine(newline);
//        }
//
/** \brief
Record information about cohorts reaching maturity in the model
*/
@param latIndex The latitude index of the grid cell in which the cohort was spawned 
@param lonIndex The longitude index of the grid cell in which the cohort was spawned 
@param timestep The model timestep in which the spawning happened 
@param birthTimestep The timestep in which the cohort was born 
@param juvenileMass The mass at which the cohort was born 
@param adultMass The maturity mass of the cohort 
@param functionalGroup The functional group of the cohort 
//         void TrackMaturity(unsigned latIndex, unsigned lonIndex, unsigned timestep, unsigned birthTimestep, double juvenileMass, double adultMass, int functionalGroup)
//        {
//            //Record data on this cohort reaching maturity
//
//            double[] MaturityRecords = new double[4];
//            MaturityRecords[0] = (double)birthTimestep;
//            MaturityRecords[1] = juvenileMass;
//            MaturityRecords[2] = adultMass;
//            MaturityRecords[3] = (double)functionalGroup;
//
//            //_Maturity[latIndex, lonIndex,timestep].Add(MaturityRecords);
//
//            // Write the time step and the abundance of the new cohort to the output file for diagnostic purposes
//            string newline = Convert.ToString(latIndex) +'\t'+ Convert.ToString(lonIndex)+'\t'+
//                Convert.ToString(timestep) + '\t' + Convert.ToString(birthTimestep) + '\t' + Convert.ToString(juvenileMass) + '\t'+
//                Convert.ToString(adultMass) + '\t' + Convert.ToString(functionalGroup);
//            SyncMaturityWriter.WriteLine(newline);
//        }
//
/** \brief
Close the output streams for the reproduction tracker
*/
//         void CloseStreams()
//        {
//            SyncMaturityWriter.Close();
//            MaturityWriter.Close();
//            SyncNewCohortWriter.Close();
//            NewCohortWriter.Close();
//        }
//
//    }
//}
#endif
