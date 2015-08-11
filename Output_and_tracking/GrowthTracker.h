#ifndef GROWTHTRACKER_H
#define GROWTHTRACKER_H
/** \file GrowthTracker.h
 * \brief the GrowthTracker header file
 */






//
//namespace Madingley
//{
/** \brief
//    /// Track the growth of a cohort in a time step
//    /// </summary>
//     class GrowthTracker
//    {
/** \brief
File to write data on growth to
*/
//        string GrowthFilename;
//
/** \brief
A streamwriter for writing out data on growth
*/
//        private StreamWriter GrowthWriter;
//        private TextWriter SyncGrowthWriter;
//
/** \brief
Set up the tracker for outputing the growth of cohorts each time step
*/
@param numTimeSteps The total number of timesteps for this simulation 
@param numLats The number of latitudes in the model grid 
@param numLons The number of longitudes in the model grid 
@param cellIndices List of indices of active cells in the model grid 
@param growthFilename The name of the file to write information about growth to 
@param outputFilesSuffix The suffix to apply to output files from this simulation 
@param outputPath The file path to write all outputs to 
@param cellIndex The index of the current cell in the list of all cells in this simulation 
//         GrowthTracker(unsigned numTimeSteps, unsigned numLats, unsigned numLons, List<unsigned[]> cellIndices, string growthFilename,
//            string outputFilesSuffix, string outputPath, int cellIndex)
//        {
//            GrowthFilename = growthFilename;
//
//            // Initialise streamwriter to output growth data
//            GrowthWriter = new StreamWriter(outputPath + GrowthFilename + outputFilesSuffix + "_Cell" + cellIndex + ".txt");
//            SyncGrowthWriter = TextWriter.Synchronized(GrowthWriter);
//            SyncGrowthWriter.WriteLine("Latitude\tLongitude\ttime_step\tCurrent_body_mass_g\tfunctional_group\tgrowth_g\tmetabolism_g\tpredation_g\therbivory_g");
//
//        }
//
/** \brief
Record the growth of the individuals in a cohort in the current time step
*/
@param latIndex The latitudinal index of the current grid cell 
@param lonIndex The longitudinal index of the current grid cell 
@param timeStep The current time step 
@param currentBodyMass The current body mass of individuals in the cohort 
@param functionalGroup The index of the functional group that the cohort belongs to 
@param netGrowth The net growth of individuals in the cohort this time step 
@param metabolism The biomass lost by individuals in this cohort through metabolism 
@param predation The biomass gained by individuals in this cohort through predation 
@param herbivory The biomass gained by individuals in this cohort through herbivory 
//         void RecordGrowth(unsigned latIndex, unsigned lonIndex, unsigned timeStep, double currentBodyMass, int functionalGroup, 
//            double netGrowth, double metabolism, double predation, double herbivory)
//        {
//            SyncGrowthWriter.WriteLine(Convert.ToString(latIndex) + '\t' + Convert.ToString(lonIndex) + '\t' + Convert.ToString(timeStep) +
//                '\t' + Convert.ToString(currentBodyMass) + '\t' + Convert.ToString(functionalGroup) + '\t' + Convert.ToString(netGrowth)+ '\t' + Convert.ToString(metabolism)+ '\t' + Convert.ToString(predation)+ '\t' + Convert.ToString(herbivory));
//        }
//
/** \brief
Closes streams for writing growth data
*/
//         void CloseStreams()
//        {
//
//            SyncGrowthWriter.Dispose();
//            GrowthWriter.Dispose();
//        }
//
//    }
//}
#endif
