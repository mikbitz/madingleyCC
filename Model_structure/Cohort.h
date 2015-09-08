#ifndef COHORT_H
#define COHORT_H
#include <Properties.h>
#include <limits.h>
/** \file Cohort.h
 * \brief the Cohort header file
 */


//namespace Madingley
//{
/** \brief Class to hold properties of a single cohort */
class Cohort
    {
    public:
//        
/** \brief Time step when the cohort was generated */
       unsigned BirthTimeStep;
/** \brief The time step at which this cohort reached maturity */
         UnsignedProperty MaturityTimeStep;
/** \brief A list of all cohort IDs ever associated with individuals in this current cohort */
        vector<long> CohortID;
/** \brief The mean juvenile mass of individuals in this cohort */
         double JuvenileMass;
/** \brief The mean mature adult mass of individuals in this cohort */
           double AdultMass;
/** \brief The mean body mass of an individual in this cohort */
         DoubleProperty IndividualBodyMass;
/** \brief Individual biomass assigned to reproductive potential */
          double IndividualReproductivePotentialMass;
/** \brief The maximum mean body mass ever achieved by individuals in this cohort */
        double MaximumAchievedBodyMass;  
/** \brief The number of individuals in the cohort */
          double CohortAbundance;
/** \brief The index of the functional group that the cohort belongs to */
//here the original code said "byte" - really necessary? not sure if unsigned char will work
         unsigned char FunctionalGroupIndex;
/** \brief Whether this cohort has ever been merged with another cohort */
        bool Merged;
/** \brief The proportion of the timestep for which this cohort is active */
         double ProportionTimeActive;
/** \brief The optimal prey body size for individuals in this cohort */
        double LogOptimalPreyBodySizeRatio;
/** \brief
Constructor for the Cohort class: assigns cohort starting properties

@param functionalGroupIndex The functional group index of the cohort being generated 
@param juvenileBodyMass The mean juvenile body mass of individuals in the cohort 
@param adultBodyMass The mean mature adult body mass of individuals in the cohort 
@param initialBodyMass The intial mean body mass of individuals in this cohort 
@param initialAbundance The intial number of individuals in this cohort 
@param optimalPreyBodySizeRatio The optimal prey body mass (as a percentage of this cohorts mass) for individuals in this cohort 
@param birthTimeStep The birth time step for this cohort 
@param nextCohortID The unique ID to assign to the next cohort created 
@param tracking Whether the process tracker is enabled 
*/
Cohort(unsigned char functionalGroupIndex, double juvenileBodyMass, double adultBodyMass, double initialBodyMass, double initialAbundance, double optimalPreyBodySizeRatio, unsigned short birthTimeStep, double proportionTimeActive, long long &nextCohortID, bool tracking)
       {
           FunctionalGroupIndex = functionalGroupIndex;
           JuvenileMass = juvenileBodyMass;
           AdultMass = adultBodyMass;
           IndividualBodyMass = initialBodyMass;
           CohortAbundance = initialAbundance;
           BirthTimeStep = birthTimeStep;
           MaturityTimeStep = std::numeric_limits<unsigned>::max();
           LogOptimalPreyBodySizeRatio = log(optimalPreyBodySizeRatio);
           MaximumAchievedBodyMass = juvenileBodyMass;
           Merged = false;
           ProportionTimeActive = proportionTimeActive;
           //if(tracking)_CohortID.Add(Convert.ToUInt32(nextCohortID));
           nextCohortID++;
        }
    };
//}
#endif
