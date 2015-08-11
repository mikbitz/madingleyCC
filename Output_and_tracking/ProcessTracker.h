#ifndef PROCESSTRACKER_H
#define PROCESSTRACKER_H
/** \file ProcessTracker.h
 * \brief the ProcessTracker header file
 */


//
//namespace Madingley
//{
/** \brief  Tracks diagnostics about the ecological processes */
     class ProcessTracker
    {
    public:
/** \brief Whether to track processes*/
//        private Boolean _TrackProcesses;
/** \brief Get or set whether to track processes*/
         bool TrackProcesses;
//        {
//            get { return _TrackProcesses; }
//            set { _TrackProcesses = value; }
//        }
//        
/** \brief Instance of the reproduction tracker within the process tracker*/
//        private ReproductionTracker  _TrackReproduction;
/** \brief Get and set the reproduction tracker*/
//         ReproductionTracker  TrackReproduction
//        {
//            get { return _TrackReproduction; }
//            set { _TrackReproduction = value; }
//        }
//
/** \brief Instance of predation tracker*/
//        private PredationTracker  _TrackPredation;
/** \brief Get and set the predation tracker*/
//         PredationTracker  TrackPredation
//        {
//            get { return _TrackPredation; }
//            set { _TrackPredation = value; }
//        }
//
/** \brief Instance of the eating tracker*/
//        private EatingTracker _TrackEating;
/** \brief Get and set the eating tracker*/
//         EatingTracker TrackEating
//        {
//            get { return _TrackEating; }
//            set { _TrackEating = value; }
//        }
//
/** \brief Instance of the growth tracker*/
//        private GrowthTracker _TrackGrowth;
/** \brief Get and set the growth tracker*/
//         GrowthTracker TrackGrowth
//        {
//            get { return _TrackGrowth; }
//            set { _TrackGrowth = value; }
//        }
//
/** \brief Instance of the mortality tracker*/
//        private MortalityTracker _TrackMortality;
/** \brief Get and set the mortality tracker*/
//         MortalityTracker TrackMortality
//        {
//            get { return _TrackMortality; }
//            set { _TrackMortality = value; }
//        }
//
/** \brief An instance of the extinction tracker*/
//        private ExtinctionTracker _TrackExtinction;
/** \brief Get and set the instance of the extinction tracker*/
//         ExtinctionTracker TrackExtinction
//        {
//            get { return _TrackExtinction; }
//            set { _TrackExtinction = value; }
//        }
//
/** \brief An instance of the metabolism tracker*/
//        private MetabolismTracker _TrackMetabolism;
/** \brief Get and set the instance of the metabolism tracker*/
//         MetabolismTracker TrackMetabolism
//        {
//            get { return _TrackMetabolism; }
//            set { _TrackMetabolism = value; }
//        }
//       
//
//
/** \brief Constructor for process tracker: Initialises the trackers for individual processes
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
@param specificLocations Whether the model is being run for specific locations  */
//         ProcessTracker(unsigned numTimesteps,
//            float[] lats, float[] lons, 
//            List<unsigned[]> cellIndices,
//            map<string,string> Filenames, 
//            Boolean trackProcesses, 
//            FunctionalGroupDefinitions cohortDefinitions, 
//            double missingValue,
//            string outputFileSuffix,
//            string outputPath, MassBinsHandler trackerMassBins,
//            Boolean specificLocations,
//            int cellIndex,
//            MadingleyModelInitialisation initialisation,
//            bool marineCell,
//            float latCellSize,
//            float lonCellSize)
//        {
//            // Initialise trackers for ecological processes
//            _TrackProcesses = trackProcesses;
//
//            if (_TrackProcesses)
//            {
//                _TrackReproduction = new ReproductionTracker(numTimesteps, (unsigned)lats.Length, (unsigned)lons.Length, cellIndices, Filenames["NewCohortsOutput"], Filenames["MaturityOutput"], outputFileSuffix, outputPath, cellIndex);
//                _TrackEating = new EatingTracker((unsigned)lats.Length, (unsigned)lons.Length, Filenames["TrophicFlowsOutput"], outputFileSuffix, outputPath, cellIndex, initialisation, marineCell);
//                _TrackGrowth = new GrowthTracker(numTimesteps, (unsigned)lats.Length, (unsigned)lons.Length, cellIndices, Filenames["GrowthOutput"], outputFileSuffix, outputPath, cellIndex);
//                _TrackMortality = new MortalityTracker(numTimesteps, (unsigned)lats.Length, (unsigned)lons.Length, cellIndices, Filenames["MortalityOutput"], outputFileSuffix, outputPath, cellIndex);
//                _TrackExtinction = new ExtinctionTracker(Filenames["ExtinctionOutput"], outputPath, outputFileSuffix, cellIndex);
//                _TrackMetabolism = new MetabolismTracker(Filenames["MetabolismOutput"], outputPath, outputFileSuffix, cellIndex);
//
//                // Initialise the predation and herbivory trackers only for runs with specific locations
//                if (specificLocations == true)
//                {
//                    _TrackPredation = new PredationTracker( numTimesteps, cellIndices, Filenames["PredationFlowsOutput"], cohortDefinitions,
//                        missingValue, outputFileSuffix, outputPath, trackerMassBins, cellIndex);
//                }
//            }
//        }
//
/** \brief Record a new cohort in the reproduction tracker
@param latIndex The latitudinal index of the current grid cell 
@param lonIndex The longitudinal index of the current grid cell 
@param timestep The current model time step 
@param offspringCohortAbundance The number of individuals in the new cohort 
@param parentCohortAdultMass The adult body mass of the parent cohort 
@param functionalGroup The functional group that the parent and offspring cohorts belong to 
@param parentCohortIDs All cohort IDs associated with the acting parent cohort 
@param offspringCohortID The cohort ID that has been assigned to the produced offspring cohort  */
//         void RecordNewCohort(unsigned latIndex, unsigned lonIndex, unsigned timestep, double offspringCohortAbundance, 
//            double parentCohortAdultMass, int functionalGroup, List<unsigned> parentCohortIDs, unsigned offspringCohortID)
//        {
//            _TrackReproduction.RecordNewCohort(latIndex, lonIndex, timestep, offspringCohortAbundance, parentCohortAdultMass, 
//                functionalGroup,parentCohortIDs,offspringCohortID);
//        }
//
/** \brief Track the maturity of a cohort in the reproduction tracker
@param latIndex The latitudinal index of the current grid cell 
@param lonIndex The longitudinal index of the current grid cell 
@param timestep The current model time step 
@param birthTimestep The birth time step of the cohort reaching maturity 
@param juvenileMass The juvenile mass of the cohort reaching maturity 
@param adultMass The adult mass of the cohort reaching maturity 
@param functionalGroup The functional group of the cohort reaching maturity  */
//         void TrackMaturity(unsigned latIndex, unsigned lonIndex, unsigned timestep, unsigned birthTimestep, double juvenileMass, double adultMass, int functionalGroup)
//        {
//            _TrackReproduction.TrackMaturity(latIndex,lonIndex,timestep,birthTimestep,juvenileMass,adultMass,functionalGroup);
//        }
//
/** \brief Track the flow of mass between trophic levels during a predation event
@param latIndex The latitudinal index of the current grid cell 
@param lonIndex The longitudinal index of the current grid cell 
@param fromFunctionalGroup The index of the functional group being eaten 
@param toFunctionalGroup The index of the functional group that the predator belongs to 
@param cohortFunctionalGroupDefinitions The functional group definitions of cohorts in the model 
@param massEaten The mass eaten during the predation event  */
//         void TrackPredationTrophicFlow(unsigned latIndex, unsigned lonIndex, int fromFunctionalGroup, int toFunctionalGroup,FunctionalGroupDefinitions cohortFunctionalGroupDefinitions, double massEaten, double predatorBodyMass, double preyBodyMass, MadingleyModelInitialisation initialisation, Boolean marineCell)
//        {
//            _TrackEating.RecordPredationTrophicFlow(latIndex, lonIndex, fromFunctionalGroup, toFunctionalGroup, cohortFunctionalGroupDefinitions, massEaten, predatorBodyMass, preyBodyMass, initialisation, marineCell);
//        }
//
/** \brief Track the flow of mass between trophic levels during a herbivory event
@param latIndex The latitudinal index of the current grid cell 
@param lonIndex The longitudinal index of the current grid cell 
@param toFunctionalGroup The index of the functional group that the predator belongs to 
@param cohortFunctionalGroupDefinitions The functional group definitions of cohorts in the model 
@param massEaten The mass eaten during the herbivory event  */
//         void TrackHerbivoryTrophicFlow(unsigned latIndex, unsigned lonIndex, int toFunctionalGroup, FunctionalGroupDefinitions cohortFunctionalGroupDefinitions, double massEaten, double predatorBodyMass, MadingleyModelInitialisation initialisation, Boolean marineCell)
//        {
//            _TrackEating.RecordHerbivoryTrophicFlow(latIndex, lonIndex, toFunctionalGroup, cohortFunctionalGroupDefinitions, massEaten, predatorBodyMass, initialisation, marineCell);
//        }
//
/** \brief Track the flow of mass between trophic levels during primary production of autotrophs
@param latIndex The latitudinal index of the current grid cell 
@param lonIndex The longitudinal index of the current grid cell 
@param massEaten The mass gained through primary production  */
//         void TrackPrimaryProductionTrophicFlow(unsigned latIndex, unsigned lonIndex, double massEaten)
//        {
//            _TrackEating.RecordPrimaryProductionTrophicFlow(latIndex, lonIndex, massEaten);
//        }
//
/** \brief Write the 
@param currentTimeStep  
@param numLats  
@param numLons   */
//         void WriteTimeStepTrophicFlows(unsigned currentTimeStep,unsigned numLats,unsigned numLons, MadingleyModelInitialisation initialisation, Boolean marineCell)
//        {
//            _TrackEating.WriteTrophicFlows(currentTimeStep, numLats, numLons, initialisation, marineCell);
//        }
//
/** \brief Track growth of individuals in a cohort using the growth tracker
@param latIndex The latitudinal index of the current grid cell 
@param lonIndex The longitudinal index of the current grid cell 
@param timeStep The current model time step 
@param currentBodyMass The current body mass of individuals in the cohort 
@param functionalGroup The funcitonal group of the cohort being tracked 
@param netGrowth The net growth of individuals in the cohort this time step 
@param metabolism The mass lost to indivduals in the cohort through metabolism 
@param predation The mass gained by individuals in the cohort through predation 
@param herbivory The mass gained by individuals in the cohort through herbivory  */
//         void TrackTimestepGrowth(unsigned latIndex, unsigned lonIndex, unsigned timeStep, double currentBodyMass, 
//            int functionalGroup, double netGrowth, double metabolism, double predation, double herbivory)
//        {
//            _TrackGrowth.RecordGrowth(latIndex, lonIndex, timeStep, currentBodyMass, functionalGroup, netGrowth, metabolism, predation,herbivory);
//        }
//
/** \brief Records the flow of mass between a prey and its predator during a predation event
@param currentTimeStep The current model time step 
@param preyBodyMass The individual body mass of the prey 
@param predatorBodyMass The individual body mass of the predator 
@param massFlow The flow of mass between predator and prey  */
//         void RecordPredationMassFlow(unsigned currentTimeStep, double preyBodyMass, double predatorBodyMass, double massFlow)
//        {
//            _TrackPredation.RecordFlow(currentTimeStep, preyBodyMass, predatorBodyMass, massFlow);
//        }
//
/** \brief Adds the mass flows from predation in the current time step to the output file and then resets the mass flow tracker
@param currentTimeStep The current model time step  */
//         void EndTimeStepPredationTracking(unsigned currentTimeStep)
//        {
//            _TrackPredation.AddTimestepFlows((int)currentTimeStep);
//            _TrackPredation.ResetPredationTracker();
//        }
//
/** \brief Records the flow of mass between primary producers and herbivores during a herbivory event
@param currentTimeStep The current model time step 
@param herbivoreBodyMass The individual body mass of the herbivore 
@param massFlow The flow of mass between the primary producer and the herbivore  */
//         void RecordHerbivoryMassFlow(unsigned currentTimeStep, double herbivoreBodyMass, double massFlow)
//        {
//            //_TrackHerbivory.RecordFlow(currentTimeStep, herbivoreBodyMass, massFlow);
//        }
//
/** \brief Adds the mass flows from herbivory in the current time step to the output file and then resets the mass flow tracker
@param currentTimeStep   */
//         void EndTimeStepHerbvioryTracking(unsigned currentTimeStep)
//        {
//            //_TrackHerbivory.AddTimestepFlows((int)currentTimeStep);
//            //_TrackHerbivory.ResetHerbivoryTracker();
//        }
//
/** \brief Record an instance of mortality in the output file
@param latIndex The latitudinal index of the current grid cell 
@param lonIndex The longitudinal index of the current grid cell 
@param birthTimeStep The time step in which this cohort was born 
@param timeStep The current model time step 
@param currentMass The current body mass of individuals in the cohort 
@param adultMass The adult mass of individuals in the cohort 
@param functionalGroup The functional group of the cohort suffering mortality 
@param cohortID The ID of the cohort suffering mortality 
@param numberDied The number of individuals dying in this mortality event 
@param mortalitySource The type of mortality causing the individuals to die  */
//         void RecordMortality(unsigned latIndex, unsigned lonIndex, unsigned birthTimeStep, unsigned timeStep, double currentMass, double adultMass, unsigned functionalGroup, unsigned cohortID, 
//            double numberDied,string mortalitySource)
//        {
//            _TrackMortality.RecordMortality(latIndex, lonIndex, birthTimeStep,
//                timeStep, currentMass, adultMass, functionalGroup, cohortID, numberDied, mortalitySource);
//        }
//
/** \brief Output the mortality profile of a cohort becoming extinct
@param cohortID The ID of the cohort becoming extinct  */
//         void OutputMortalityProfile(unsigned cohortID)
//        {
//            _TrackMortality.OutputMortalityProfile(cohortID);
//        }
//
/** \brief Record the extinction of a cohort
@param latIndex The latitudinal index of the current grid cell 
@param lonIndex The longitudinal index of the current grid cell 
@param currentTimeStep THe current time step 
@param merged Whether the cohort becoming extinct has ever been merged 
@param cohortIDs The IDs of all cohorts that have contributed individuals to the cohort going extinct  */
//         void RecordExtinction(unsigned latIndex, unsigned lonIndex,unsigned currentTimeStep,bool merged,List<unsigned>cohortIDs)
//        {
//            _TrackExtinction.RecordExtinction(latIndex, lonIndex, currentTimeStep, merged, cohortIDs);
//        }
//
/** \brief Tracks the mass lost by individuals in a cohort in a time step through metabolism
@param latIndex The latitudinal index of the current grid cell 
@param lonIndex The longitudinal index of the current grid cell 
@param timeStep The current model time step 
@param currentBodyMass The body mass of individuals in the acting cohort 
@param functionalGroup The functional group index of the acting cohort 
@param temperature The ambient temperature in the grid cell 
@param metabolicLoss The mass lost by individuals through metabolism  */
//         void TrackTimestepMetabolism(unsigned latIndex, unsigned lonIndex, unsigned timeStep, double currentBodyMass, 
//            int functionalGroup, double temperature, double metabolicLoss)
//        {
//            _TrackMetabolism.RecordMetabolism(latIndex, lonIndex, timeStep, currentBodyMass, functionalGroup, temperature, metabolicLoss);
//        }
//
//
/** \brief Close all tracker streams*/
//		 void CloseStreams(Boolean SpecificLocations)
//        {
//            _TrackReproduction.CloseStreams();
//            _TrackEating.CloseStreams();
//            _TrackGrowth.CloseStreams();
//            _TrackMetabolism.CloseStreams();
//            //_TrackNPP.CloseStreams();
//            if (SpecificLocations == true)
//            {
//                _TrackPredation.CloseStreams();
//            }
//        }
    };
//}
#endif
