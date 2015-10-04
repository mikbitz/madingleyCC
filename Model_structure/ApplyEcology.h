#ifndef APPLYECOLOGY_H
#define APPLYECOLOGY_H
#include<ProcessTracker.h>
/** \file ApplyEcology.h
 * \brief the ApplyEcology header file
 */







//
//namespace Madingley
//{
/** \brief Class for applying changes from the ecological processes to the properties of the acting cohort and to the environment */
    class ApplyEcology
   {
public:
/** \brief
Apply all updates from the ecological processes to the properties of the acting cohort and to the environment

@param gridCellCohorts The cohorts in the current grid cell 
@param actingCohort The location of the acting cohort in the jagged array of grid cell cohorts 
@param cellEnvironment The environment in the current gird cell 
@param deltas The sorted list to track changes in biomass and abundance of the acting cohort in this grid cell 
@param currentTimestep The current model time step 
@param tracker A process tracker */
void UpdateAllEcology(GridCellCohortHandler& gridCellCohorts, vector<int>& actingCohort, map<string, vector<double>>& cellEnvironment, map<string, map<string, double>>& 
           deltas, unsigned currentTimestep, ProcessTracker& tracker)
       {
          // Apply cohort abundance changes
           UpdateAbundance(gridCellCohorts, actingCohort, deltas); 
           // Apply cohort biomass changes
           UpdateBiomass(gridCellCohorts, actingCohort, deltas, currentTimestep, tracker, cellEnvironment);
           // Apply changes to the environmental biomass pools
           UpdatePools(cellEnvironment, deltas);
       }

/** \brief
Update the abundance of the acting cohort according to the delta abundances from the ecological processes

@param gridCellCohorts The cohorts in the current grid cell 
@param actingCohort The location of the acting cohort in the jagged array of grid cell cohorts 
@param deltas The sorted list to track changes in biomass and abundance of the acting cohort in this grid cell */
void UpdateAbundance(GridCellCohortHandler& gridCellCohorts, vector<int>& actingCohort, map<string,map<string, double>>& deltas)
       {
           // Extract the abundance deltas from the sorted list of all deltas
           map<string, double> deltaAbundance = deltas["abundance"];
                     
           // Variable to calculate net abundance change to check that cohort abundance will not become negative
           double NetAbundanceChange = 0.0;

           // Loop over all abundance deltas
           for (auto d : deltaAbundance )
           {
               // Update net abundance change
               NetAbundanceChange += d.second;
           }
           // Check that cohort abundance will not become negative
           assert((gridCellCohorts[actingCohort].CohortAbundance + NetAbundanceChange) >= 0 && "Cohort abundance < 0");

           //Loop over all keys in the abundance deltas sorted list
           for (auto d : deltaAbundance)
           {
               // Update the abundance of the acting cohort
               gridCellCohorts[actingCohort].CohortAbundance += d.second;
               // Reset the current delta abundance to zero
               d.second = 0.0;
           }

       }

/** \brief
Update the individual and reproductive body masses of the acting cohort according to the delta biomasses from the ecological processes

@param gridCellCohorts The cohorts in the current grid cell 
@param actingCohort The position of the acting cohort in the jagged array of grid cell cohorts 
@param deltas The sorted list to track changes in biomass and abundance of the acting cohort in this grid cell 
@param currentTimestep The current model time step 
@param tracker A process tracker 
@param cellEnvironment The cell environment */
void UpdateBiomass(GridCellCohortHandler& gridCellCohorts, vector<int>& actingCohort, map<string,map<string, double>>& deltas, 
           unsigned currentTimestep, ProcessTracker& tracker, map<string,vector<double>>& cellEnvironment)
       {
           // Extract the biomass deltas from the sorted list of all deltas
           map<string, double> deltaBiomass = deltas["biomass"];

//            if (tracker.TrackProcesses)
//            {
//                // Calculate net growth of individuals in this cohort
//                double growth = deltaBiomass["predation"] + deltaBiomass["herbivory"] + deltaBiomass["metabolism"];
//                tracker.TrackTimestepGrowth((unsigned)cellEnvironment["LatIndex"][0], (unsigned)cellEnvironment["LonIndex"][0], currentTimestep,
//                    gridCellCohorts[actingCohort].IndividualBodyMass, gridCellCohorts[actingCohort].FunctionalGroupIndex, growth, deltaBiomass["metabolism"],deltaBiomass["predation"],deltaBiomass["herbivory"]);
//                  
//            }
//

           // Variable to calculate net biomass change to check that cohort individual body mass will not become negative
           double NetBiomass = 0.0;

           // Loop over all biomass deltas
           for (auto d : deltaBiomass)
           {
               // Update net biomass change
               NetBiomass += d.second;
           }
           //cout<<"oower "<<NetBiomass<<endl;
           double BiomassCheck=0.0;
           bool NetToBeApplied = true;
           // If cohort abundance is greater than zero, then check that the calculated net biomass will not make individual body mass become negative
           if (gridCellCohorts[actingCohort].CohortAbundance > 0)
           {

               BiomassCheck = gridCellCohorts[actingCohort].IndividualBodyMass() + NetBiomass;
               if (BiomassCheck < 0){
                  cout<<"Biomass going negative, acting cohort: "<< actingCohort[0] << ", "<< actingCohort[1];
                  exit(1); 
               }
           }

           //Loop over all keys in the abundance deltas sorted list
           for (auto d :deltaBiomass)
           {
               // If cohort abundance is zero, then set cohort individual body mass to zero and reset the biomass delta to zero, 
               // otherwise update cohort individual body mass and reset the biomass delta to zero
               if (gridCellCohorts[actingCohort].CohortAbundance == 0)
               {
                   gridCellCohorts[actingCohort].IndividualBodyMass = 0.0;
                   d.second = 0.0;
               }
               else
               {
                   if (NetToBeApplied)
                   {
                       gridCellCohorts[actingCohort].IndividualBodyMass = gridCellCohorts[actingCohort].IndividualBodyMass() +NetBiomass;
                       NetToBeApplied = false;
                   }

                   //gridCellCohorts[actingCohort].IndividualBodyMass += deltaBiomass[key];
                   d.second = 0.0;
               }
           }

           // Check that individual body mass is still greater than zero
           assert(gridCellCohorts[actingCohort].IndividualBodyMass() >= 0 && "biomass < 0");

           // If the current individual body mass is the largest that has been achieved by this cohort, then update the maximum achieved
           // body mass tracking variable for the cohort
           if (gridCellCohorts[actingCohort].IndividualBodyMass() > gridCellCohorts[actingCohort].MaximumAchievedBodyMass)
               gridCellCohorts[actingCohort].MaximumAchievedBodyMass = gridCellCohorts[actingCohort].IndividualBodyMass();

           // Extract the reproductive biomass deltas from the sorted list of all deltas
           map<string, double> deltaReproductiveBiomass = deltas["reproductivebiomass"];

           // Variable to calculate net reproductive biomass change to check that cohort individual body mass will not become negative
           double NetReproductiveBiomass = 0.0;

           // Loop over all reproductive biomass deltas
           for (auto d : deltaReproductiveBiomass)
           {
               // Update net reproductive biomass change
               NetReproductiveBiomass += d.second;
           }

           //Loop over all keys in the abundance deltas sorted list
           for (auto d : deltaReproductiveBiomass)
           {
               // If cohort abundance is zero, then set cohort reproductive body mass to zero and reset the biomass delta to zero, 
               // otherwise update cohort reproductive body mass and reset the biomass delta to zero
               if (gridCellCohorts[actingCohort].CohortAbundance == 0)
               {
                   gridCellCohorts[actingCohort].IndividualReproductivePotentialMass = 0.0;
                   d.second = 0.0;
               }
               else
               {
                   gridCellCohorts[actingCohort].IndividualReproductivePotentialMass += d.second;
                   d.second = 0.0;
               }
           }
           
           //Note that maturity time step is set in TReproductionBasic

       }



/** \brief
Update the organic and respiratory biomass pools according to the relevant deltas from the ecological processes

@param cellEnvironment The environment of the current gird cell 
@param deltas The sorted list to track changes in biomass and abundance of the acting cohort in this grid cell */
void UpdatePools(map<string, vector<double>>& cellEnvironment, map<string,map<string,double>>& deltas)
       {
//            // Extract the organic pool deltas from the sorted list of all deltas
//            map<string, double> DeltaOrganicPool = deltas["organicpool"];
//
//            // Get all the keys from the oragnic pool deltas sorted list
//            string[] KeyStrings = DeltaOrganicPool.Keys.ToArray();
//
//            // Loop over all keys in the organic pool deltas sorted list
//            foreach (var key in KeyStrings)
//            {
//                // Check that the delta value is not negative
//                //Debug.Assert(DeltaOrganicPool[key] >= 0.0, "A delta value for the organic pool is negative");
//                // Update the organic pool biomass
//                cellEnvironment["Organic Pool"][0] += DeltaOrganicPool[key];
//                //Reset the delta value to zero
//                DeltaOrganicPool[key] = 0.0;
//
//            }
//
//            // Extract the respiratory pool deltas from the sorted list of all deltas
//            map<string, double> DeltaRespiratoryPool = deltas["respiratoryCO2pool"];
//
//            // Get all the keys from the respiratory pool deltas sorted list
//            KeyStrings = DeltaRespiratoryPool.Keys.ToArray();
//
//            // Loop over all keys in the respiratory pool deltas sorted list
//            foreach (var key in KeyStrings)
//            {
//                // Check that the delta value is not negative
//                Debug.Assert(DeltaRespiratoryPool[key] >= 0.0, "A delta value for the respiratory CO2 pool is negative");
//                // Update the respiratory CO2 pool
//                cellEnvironment["Respiratory CO2 Pool"][0] += DeltaRespiratoryPool[key];
//                // Reset the delta value to zero
//                DeltaRespiratoryPool[key] = 0.0;
//
//            }
       }

   };
//}
#endif
