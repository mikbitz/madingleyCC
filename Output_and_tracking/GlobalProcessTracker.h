#ifndef GLOBALPROCESSTRACKER_H
#define GLOBALPROCESSTRACKER_H
/** \file GlobalProcessTracker.h
 * \brief the GlobalProcessTracker header file
 */





//
//namespace Madingley
//{
//     class GlobalProcessTracker
//    {
//                /// <summary>
Whether to track processes
*/
//        private Boolean _TrackProcesses;
/** \brief
Get or set whether to track processes
*/
//         Boolean TrackProcesses
//        {
//            get { return _TrackProcesses; }
//            set { _TrackProcesses = value; }
//        }
//
//        //An instance of the NPP tracker
//        private GlobalNPPTracker _TrackNPP;
//
//        //Get and set this instance of the NPP Tracker
//         GlobalNPPTracker TrackNPP
//        {
//            get { return _TrackNPP; }
//            set { _TrackNPP = value; }
//        }
//
/** \brief
Constructor for process tracker: Initialises the trackers for individual processes
*/
@param numTimesteps The number of time steps in the model 
@param lats The latitudes of active grid cells in the model 
@param lons The longitudes of active grid cells in the model 
@param cellIndices List of indices of active cells in the model grid 
@param Filenames The filenames of the output files to write the tracking results to 
@param trackProcesses Whether to track processes 
@param cohortDefinitions The definitions for cohort functional groups in the model 
@param missingValue The missing value to use in process tracking output files 
@param outputFileSuffix The suffix to be applied to output files from process tracking 
@param outputPath The path to the folder to be used for process tracking outputs 
@param trackerMassBins The mass bins to use for categorising output data in the process trackers 
@param specificLocations Whether the model is being run for specific locations 
//         GlobalProcessTracker(unsigned numTimesteps,
//            float[] lats, float[] lons,
//            List<unsigned[]> cellIndices,
//            map<string, string> Filenames,
//            Boolean trackProcesses,
//            FunctionalGroupDefinitions cohortDefinitions,
//            double missingValue,
//            string outputFileSuffix,
//            string outputPath, MassBinsHandler trackerMassBins,
//            Boolean specificLocations,
//            MadingleyModelInitialisation initialisation,
//            float latCellSize,
//            float lonCellSize)
//        {
//            // Initialise trackers for ecological processes
//            _TrackProcesses = trackProcesses;
//
//            if (_TrackProcesses)
//            {
//                _TrackNPP = new GlobalNPPTracker(outputPath, lats.Length, lons.Length, lats, lons, latCellSize, lonCellSize, (int)numTimesteps);
//
//            }
//        }
//
//         void RecordNPP(unsigned latIndex, unsigned lonIndex, double val)
//        {
//            _TrackNPP.RecordNPPValue(latIndex, lonIndex, val);
//        }
//
//         void StoreNPPGrid(unsigned t)
//        {
//            _TrackNPP.StoreNPPGrid(t);
//        }
//
//         void CloseNPPFile()
//        {
//            _TrackNPP.CloseNPPFile();
//
//        }
//
//    }
//}
#endif
