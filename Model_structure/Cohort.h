#ifndef COHORT_H
#define COHORT_H
#include <Properties.h>
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
/** \brief
Time step when the cohort was generated
*/
       unsigned BirthTimeStep;
/** \brief
Get time step when the cohort was generated
*/
//         unsigned BirthTimeStep { get { return _BirthTimeStep; } }
//
/** \brief
The time step at which this cohort reached maturity
*/
//        private unsigned _MaturityTimeStep;
/** \brief
Get and set the time step at which this cohort reached maturity
*/
         UnsignedProperty MaturityTimeStep;// { get { return _MaturityTimeStep; }    set { _MaturityTimeStep = value; }}
//
/** \brief
A list of all cohort IDs ever associated with individuals in this current cohort
*/
//        private List<UInt32> _CohortID= new List<UInt32>();
/** \brief
Get the list of all cohort IDs ever associated with individuals in this current cohort
*/
//         List<UInt32> CohortID {  get { return _CohortID; } }
//      
//
/** \brief
The mean juvenile mass of individuals in this cohort
*/
         double JuvenileMass;
/** \brief
Get the mean juvenile mass of individuals in this cohort
*/
//         double JuvenileMass { get { return _JuvenileMass; } }
//
/** \brief
The mean mature adult mass of individuals in this cohort
*/
           double AdultMass;
/** \brief
Get the mean mature adult mass of individuals in this cohort
*/
//         double AdultMass { get { return _AdultMass; } }
//
/** \brief
The mean body mass of an individual in this cohort
*/
//        private double _IndividualBodyMass;
/** \brief
Get or set the mean body mass of an individual in this cohort
*/
         DoubleProperty IndividualBodyMass;
//        {
//            get { return _IndividualBodyMass; }
//            set { _IndividualBodyMass = value; }
//        }
//
/** \brief
Individual biomass assigned to reproductive potential
*/
          double IndividualReproductivePotentialMass;
/** \brief
Get or set the individual biomass assigned to reproductive potential
*/
//         double IndividualReproductivePotentialMass
//        {
//            get { return _IndividualReproductivePotentialMass; }
//            set { _IndividualReproductivePotentialMass = value; }
//        }
//
/** \brief
The maximum mean body mass ever achieved by individuals in this cohort
*/
        double MaximumAchievedBodyMass;
/** \brief
Get or set the maximum mean body mass ever achieved by individuals in this cohort
*/
//         double MaximumAchievedBodyMass
//        {
//            get { return _MaximumAchievedBodyMass; }
//            set { _MaximumAchievedBodyMass = value; }
//        }
//        
/** \brief
The number of individuals in the cohort
*/
          double CohortAbundance;
/** \brief
Get or set the number of individuals in the cohort
*/
//         double CohortAbundance
//        {
//            get { return _CohortAbundance; }
//            set { _CohortAbundance = value; }
//        }
//
/** \brief
The index of the functional group that the cohort belongs to
*/
//        private byte _FunctionalGroupIndex;
/** \brief
Get the index of the functional group that the cohort belongs to
*/
//here the original code said "byte" - really necessary? not sure if unsigned char will work
         unsigned char FunctionalGroupIndex;// { get { return _FunctionalGroupIndex; } }
//
/** \brief
Whether this cohort has ever been merged with another cohort
*/
//        private Boolean _Merged;
/** \brief
Get or set whether this cohort has ever been merged with another cohort
*/
//         Boolean Merged
//        {
//            get { return _Merged; }
//            set { _Merged = value; }
//        }
//
/** \brief
The proportion of the timestep for which this cohort is active
*/
//        private double _ProportionTimeActive;
/** \brief
Get and set the proportion of time for which this cohort is active
*/
         double ProportionTimeActive;
//        {
//            get { return _ProportionTimeActive; }
//            set { _ProportionTimeActive = value; }
//        }
//        
//
/** \brief
The optimal prey body size for individuals in this cohort
*/
        double LogOptimalPreyBodySizeRatio;
/** \brief
Get and set the optimal prey body size for individuals in this cohort
*/
//         double LogOptimalPreyBodySizeRatio 
//        {
//            get { return _LogOptimalPreyBodySizeRatio ; }
//            set { _LogOptimalPreyBodySizeRatio = value; }
//        }
//        
//
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
//            _FunctionalGroupIndex = functionalGroupIndex;
//            _JuvenileMass = juvenileBodyMass;
//            _AdultMass = adultBodyMass;
//            _IndividualBodyMass = initialBodyMass;
//            _CohortAbundance = initialAbundance;
//            _BirthTimeStep = birthTimeStep;
//            _MaturityTimeStep = unsigned.MaxValue;
//            _LogOptimalPreyBodySizeRatio = Math.Log(optimalPreyBodySizeRatio);
//            _MaximumAchievedBodyMass = juvenileBodyMass;
//            _Merged = false;
//            _ProportionTimeActive = proportionTimeActive;
//            if(tracking)_CohortID.Add(Convert.ToUInt32(nextCohortID));
//            nextCohortID++;
        }
    };
//}
#endif
