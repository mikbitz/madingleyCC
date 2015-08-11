#ifndef DISPERSAL_H
#define DISPERSAL_H
#include <string>
#include <IEcologicalProcessAcrossGridCells.h>
#include <IDispersalImplementation.h>
#include <MadingleyModelInitialisation.h>
#include <GridCellCohortHandler.h>
#include <limits>
#include <TAdvectiveDispersal.h>
#include <TResponsiveDispersal.h>
#include <TDiffusiveDispersal.h>
using namespace std;
/** \file Dispersal.h
 * \brief the Dispersal header file
 */


//
//namespace Madingley
//{
/** \brief Performs dispersal */
class Dispersal : IEcologicalProcessAcrossGridCells
    {
/** \brief The available implementations of the dispersal process */
          map<string, IDispersalImplementation*> Implementations;
//
/** \brief Threshold (g) below which a marine individual is considered to be planktonic (i.e. cannot swim against the currents). Currently set to 10mg. */
          double PlanktonThreshold;
//
/** \brief Constructor for Dispersal: fills the list of available implementations of dispersal */
    public:
    Dispersal(bool DrawRandomly, string globalModelTimeStepUnit, MadingleyModelInitialisation modelInitialisation)
        {
//Class member does not need initialising
            // Initialise the list of dispersal implementations
//            Implementations = new map<string, IDispersalImplementation>();
//
           // Add the basic advective dispersal implementation to the list of implementations
           AdvectiveDispersal* AdvectiveDispersalImplementation = new AdvectiveDispersal(globalModelTimeStepUnit, DrawRandomly);
           Implementations["basic advective dispersal"]= AdvectiveDispersalImplementation;

           // Add the basic advective dispersal implementation to the list of implementations
           DiffusiveDispersal* DiffusiveDispersalImplementation = new DiffusiveDispersal(globalModelTimeStepUnit, DrawRandomly);
           Implementations["basic diffusive dispersal"]=DiffusiveDispersalImplementation;

           // Add the basic advective dispersal implementation to the list of implementations
           ResponsiveDispersal* ResponsiveDispersalImplementation = new ResponsiveDispersal(globalModelTimeStepUnit, DrawRandomly);
           Implementations["basic responsive dispersal"]= ResponsiveDispersalImplementation;

           // Get the weight threshold below which organisms are dispersed planktonically
           PlanktonThreshold = modelInitialisation.PlanktonDispersalThreshold();
        }

/** \brief Run dispersal */
         void RunCrossGridCellEcologicalProcess(vector<unsigned> cellIndex, ModelGrid gridForDispersal, bool dispersalOnly, FunctionalGroupDefinitions madingleyCohortDefinitions, FunctionalGroupDefinitions madingleyStockDefinitions, unsigned currentMonth)
        {
/** \brief A cohort handler to temporarily hold the cohorts in each grid cell */
        GridCellCohortHandler WorkingGridCellCohorts;
//
//            // Get the lat and lon indices
            unsigned ii = cellIndex[0];
            unsigned jj = cellIndex[1];
//
//            // A boolean to check that the environmental layer exists
            bool varExists;
//            
//            // Check to see if the cell is marine
            double CellRealm = gridForDispersal.GetEnviroLayer("Realm", 0, ii, jj, varExists);
//
//            // Go through all of the cohorts in turn and see if they disperse
            WorkingGridCellCohorts = gridForDispersal.GetGridCellCohorts(ii, jj);
//                    
//            // Loop through cohorts, and perform dispersal according to cohort type and status
            for (int kk = 0; kk < WorkingGridCellCohorts.size(); kk++)
            {
//                // Work through the list of cohorts
                for (int ll = 0; ll < WorkingGridCellCohorts[kk].size(); ll++)
                {
//                    // Check to see if the cell is marine and the cohort type is planktonic
                    if (CellRealm == 2.0 && 
                        ((madingleyCohortDefinitions.GetTraitNames("Mobility", WorkingGridCellCohorts[kk][ll].FunctionalGroupIndex) == "planktonic") || (WorkingGridCellCohorts[kk][ll].IndividualBodyMass() <= PlanktonThreshold)))                    
                    {   
//                        // Run advective dispersal
                        Implementations["basic advective dispersal"]->RunDispersal(cellIndex, gridForDispersal, WorkingGridCellCohorts[kk][ll], kk, ll, currentMonth);
//                        
                    }
//                    // Otherwise, if mature do responsive dispersal
                    else if (WorkingGridCellCohorts[kk][ll].MaturityTimeStep() < std::numeric_limits<unsigned>::max())
                    {
//                        // Run responsive dispersal
                        Implementations["basic responsive dispersal"]->RunDispersal(cellIndex, gridForDispersal, WorkingGridCellCohorts[kk][ll], kk, ll, currentMonth);
                    }
//                    // If the cohort is immature, run diffusive dispersal
                    else
                    {
                        Implementations["basic diffusive dispersal"]->RunDispersal(cellIndex, gridForDispersal, WorkingGridCellCohorts[kk][ll], kk, ll, currentMonth);
                    }
                }
            }
//      
//            // IF THE CELL IS MARINE, RUN ADVECTIVE DISPERSAL FOR THE PHYTOPLANKTON STOCK AS WELL IN v1
        }
    };
//}
#endif
