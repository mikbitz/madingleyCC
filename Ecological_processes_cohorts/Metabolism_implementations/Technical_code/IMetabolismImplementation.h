#ifndef IMETABOLISMIMPLEMENTATION_H
#define IMETABOLISMIMPLEMENTATION_H
#include <Properties.h>
#include <GridCellCohortHandler.h>
#include <GridCellStockHandler.h>
#include <FunctionalGroupDefinitions.h>
/** \file IMetabolismImplementation.h
 * \brief the IMetabolismImplementation header file
 */

//
//namespace Madingley
//{
/** \brief Interface for implementations of the ecological process of metabolism */
class IMetabolismImplementation
    {
    public:
/** \brief Time units associated with the formulation of metabolism */
        //StringProperty TimeUnitImplementation;
/** \brief Scalar to convert from time units associated with metabolism to the global model time step unit */
        //DoubleProperty DeltaT;
/** \brief Calculate the biomass lost through metabolism and update the relevant deltas for the acting cohort
@param gridCellCohorts The cohorts in the current grid cell 
@param gridCellStocks The stocks in the current grid cell 
@param actingCohort The position of the acting cohort in the jagged array of grid cell cohorts 
@param cellEnvironment The environment in the current grid cell 
@param deltas The sorted list to track changes in biomass and abundance of the acting cohort in this grid cell 
@param madingleyCohortDefinitions The definitions for cohort functional groups in the model 
@param madingleyStockDefinitions The definitions for stock functional groups in the model 
@param currentTimestep The current model time step 
@param currentMonth The current month in the model */
      virtual void RunMetabolism(GridCellCohortHandler& gridCellCohorts, GridCellStockHandler& gridCellStocks, vector<int>& actingCohort, 
           map<string, vector<double> >& cellEnvironment, map<string, map<string, double>>& deltas, 
           FunctionalGroupDefinitions& madingleyCohortDefinitions, FunctionalGroupDefinitions& madingleyStockDefinitions, 
           unsigned currentTimestep, unsigned currentMonth){;}
    };
//}
#endif
