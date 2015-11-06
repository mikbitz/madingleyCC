#ifndef APPLYECOLOGY_H
#define APPLYECOLOGY_H
#include <ProcessTracker.h>
#include <assert.h>

#include "Cohort.h"
/** \file ApplyEcology.h
 * \brief the ApplyEcology header file
 */


/** \brief Class for applying changes from the ecological processes to the properties of the acting cohort and to the environment */
class ApplyEcology {
public:
    //----------------------------------------------------------------------------------------------
    //Methods
    //----------------------------------------------------------------------------------------------

    //----------------------------------------------------------------------------------------------
    /** \brief  Apply all updates from the ecological processes to the properties of the acting cohort and to the environment
    @param gridCellCohorts The cohorts in the current grid cell 
    @param actingCohort The location of the acting cohort in the jagged array of grid cell cohorts 
    @param cellEnvironment The environment in the current gird cell 
    @param deltas The sorted list to track changes in biomass and abundance of the acting cohort in this grid cell 
    @param currentTimestep The current model time step 
    @param tracker A process tracker */
    void UpdateAllEcology(GridCellCohortHandler& gridCellCohorts, Cohort& actingCohort, map<string, vector<double>>&cellEnvironment, map<string, map<string, double>>&
            deltas, unsigned currentTimestep) {
        // Apply cohort abundance changes
        UpdateAbundance(gridCellCohorts, actingCohort, deltas);
        // Apply cohort biomass changes
        UpdateBiomass(gridCellCohorts, actingCohort, deltas, currentTimestep, cellEnvironment);
        // Apply changes to the environmental biomass pools
        UpdatePools(cellEnvironment, deltas);
    }
    //----------------------------------------------------------------------------------------------
    /** \brief  Update the abundance of the acting cohort according to the delta abundances from the ecological processes
    @param gridCellCohorts The cohorts in the current grid cell 
    @param actingCohort The location of the acting cohort in the jagged array of grid cell cohorts 
    @param deltas The sorted list to track changes in biomass and abundance of the acting cohort in this grid cell */
    void UpdateAbundance(GridCellCohortHandler& gridCellCohorts, Cohort& actingCohort, map<string, map<string, double>>&deltas) {
        // Extract the abundance deltas from the sorted list of all deltas
        map<string, double> deltaAbundance = deltas["abundance"];

        // Variable to calculate net abundance change to check that cohort abundance will not become negative
        double NetAbundanceChange = 0.0;

        // Loop over all abundance deltas
        for (auto d : deltaAbundance) {
            // Update net abundance change
            NetAbundanceChange += d.second;
        }
        // Check that cohort abundance will not become negative
        assert((actingCohort.CohortAbundance + NetAbundanceChange) >= 0 && "Cohort abundance < 0");

        //Loop over all keys in the abundance deltas sorted list
        for (auto d : deltaAbundance) {
            // Update the abundance of the acting cohort
            actingCohort.CohortAbundance += d.second;
            // Reset the current delta abundance to zero
            d.second = 0.0;
        }

    }
    //----------------------------------------------------------------------------------------------
    /** \brief  Update the individual and reproductive body masses of the acting cohort according to the delta biomasses from the ecological processes
    @param gridCellCohorts The cohorts in the current grid cell 
    @param actingCohort The position of the acting cohort in the jagged array of grid cell cohorts 
    @param deltas The sorted list to track changes in biomass and abundance of the acting cohort in this grid cell 
    @param currentTimestep The current model time step 
    @param tracker A process tracker 
    @param cellEnvironment The cell environment */
    void UpdateBiomass(GridCellCohortHandler& gridCellCohorts, Cohort& actingCohort, map<string, map<string, double>>&deltas,
            unsigned currentTimestep,  map<string, vector<double>>&cellEnvironment) {
        // Extract the biomass deltas from the sorted list of all deltas
        map<string, double> deltaBiomass = deltas["biomass"];

        // Variable to calculate net biomass change to check that cohort individual body mass will not become negative
        double NetBiomass = 0.0;

        // Loop over all biomass deltas
        for (auto d : deltaBiomass) {
            // Update net biomass change
            NetBiomass += d.second;
        }
        double BiomassCheck = 0.0;
        bool NetToBeApplied = true;
        // If cohort abundance is greater than zero, then check that the calculated net biomass will not make individual body mass become negative
        if (actingCohort.CohortAbundance > 0) {

            BiomassCheck = actingCohort.IndividualBodyMass + NetBiomass;
            if (BiomassCheck < 0) {
                cout << "Biomass going negative, acting cohort: " << actingCohort.FunctionalGroupIndex << ", " << actingCohort.ID;
                exit(1);
            }
        }

        //Loop over all keys in the abundance deltas sorted list
        for (auto d : deltaBiomass) {
            // If cohort abundance is zero, then set cohort individual body mass to zero and reset the biomass delta to zero, 
            // otherwise update cohort individual body mass and reset the biomass delta to zero
            if (actingCohort.CohortAbundance == 0) {
                actingCohort.IndividualBodyMass = 0.0;
                d.second = 0.0;
            } else {
                if (NetToBeApplied) {
                    actingCohort.IndividualBodyMass = actingCohort.IndividualBodyMass + NetBiomass;
                    NetToBeApplied = false;
                }

                //actingCohort.IndividualBodyMass += deltaBiomass[key];
                d.second = 0.0;
            }
        }

        // Check that individual body mass is still greater than zero
        assert(actingCohort.IndividualBodyMass >= 0 && "biomass < 0");

        // If the current individual body mass is the largest that has been achieved by this cohort, then update the maximum achieved
        // body mass tracking variable for the cohort
        if (actingCohort.IndividualBodyMass > actingCohort.MaximumAchievedBodyMass)
            actingCohort.MaximumAchievedBodyMass = actingCohort.IndividualBodyMass;

        // Extract the reproductive biomass deltas from the sorted list of all deltas
        map<string, double> deltaReproductiveBiomass = deltas["reproductivebiomass"];

        // Variable to calculate net reproductive biomass change to check that cohort individual body mass will not become negative
        double NetReproductiveBiomass = 0.0;

        // Loop over all reproductive biomass deltas
        for (auto d : deltaReproductiveBiomass) {
            // Update net reproductive biomass change
            NetReproductiveBiomass += d.second;
        }

        //Loop over all keys in the abundance deltas sorted list
        for (auto d : deltaReproductiveBiomass) {
            // If cohort abundance is zero, then set cohort reproductive body mass to zero and reset the biomass delta to zero, 
            // otherwise update cohort reproductive body mass and reset the biomass delta to zero
            if (actingCohort.CohortAbundance == 0) {
                actingCohort.IndividualReproductivePotentialMass = 0.0;
                d.second = 0.0;
            } else {
                actingCohort.IndividualReproductivePotentialMass += d.second;
                d.second = 0.0;
            }
        }
        //Note that maturity time step is set in TReproductionBasic
    }
    //----------------------------------------------------------------------------------------------
    /** \brief Update the organic and respiratory biomass pools according to the relevant deltas from the ecological processes

    @param cellEnvironment The environment of the current grid cell 
    @param deltas The sorted list to track changes in biomass and abundance of the acting cohort in this grid cell */
    void UpdatePools(map<string, vector<double>>&cellEnvironment, map<string, map<string, double>>&deltas) {
        
                    // Loop over all keys in the organic pool deltas sorted list
                    for (auto &D : deltas["organicpool"])
                    {
                        // Check that the delta value is not negative
                        assert(D.second >= 0.0 && "A delta value for the organic pool is negative");
                        // Update the organic pool biomass
                        cellEnvironment["Organic Pool"][0] += D.second;
                        //Reset the delta value to zero
                        D.second = 0.0;
        
                    }

                    // Loop over all keys in the respiratory pool deltas sorted list
                    for (auto &D : deltas["respiratoryCO2pool"])
                    {
                        // Check that the delta value is not negative
                        assert(D.second >= 0.0 && "A delta value for the respiratory CO2 pool is negative");
                        // Update the respiratory CO2 pool
                        cellEnvironment["Respiratory CO2 Pool"][0] += D.second;
                        // Reset the delta value to zero
                        D.second = 0.0;
        
                    }
    }
    //----------------------------------------------------------------------------------------------
};
#endif
