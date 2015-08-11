#ifndef TBACKGROUNDMORTALITY_H
#define TBACKGROUNDMORTALITY_H
/** \file TBackgroundMortality.h
 * \brief the TBackgroundMortality header file
 */

//

//
//namespace Madingley
//{
/** \brief A formulation of the process of background mortality, i.e. mortality from disease, accidents and other random events*/
class BackgroundMortality : public IMortalityImplementation
    {
/** \brief The time units associated with this implementation of dispersal*/
           const string TimeUnitImplementation = "Day";
/** \brief Cohort background mortality rate - the proportion of individuals dying in a time step */
           const double MortalityRate = 0.001;
/** \brief Scalar to convert from the time step units used by this formulation of dispersal to global model time step units */
           DoubleProperty DeltaT;
/** \brief Include Utility class */
           UtilityFunctions Utilities;
    public:
/** \brief Constructor for background mortality: assigns all parameter values*/
BackgroundMortality(string globalModelTimeStepUnit)
        {
           // Calculate the scalar to convert from the time step units used by this implementation of mortality to the global model time step units
            DeltaT = Utilities.ConvertTimeUnits(globalModelTimeStepUnit, TimeUnitImplementation);
        }
/** \brief Calculate the rate of individuals in a cohort that die from background mortality in a model time step
@param gridCellCohorts The cohorts in the current grid cell 
@param actingCohort The position of the acting cohort in the jagged array of grid cell cohorts 
@param deltas The sorted list to track changes in biomass and abundance of the acting cohort in this grid cell 
@param currentTimestep The current model time step 
@return The rate of individuals in the cohort that die from background mortality
*/
double CalculateMortalityRate(GridCellCohortHandler gridCellCohorts, vector<int> actingCohort, double bodyMassIncludingChangeThisTimeStep, map<string, map<string, double> > deltas, unsigned currentTimestep)
       {
           // Convert from mortality rate per mortality formulation time step to mortality rate per model time step
           return MortalityRate * DeltaT();
       }
   };

//}
#endif
