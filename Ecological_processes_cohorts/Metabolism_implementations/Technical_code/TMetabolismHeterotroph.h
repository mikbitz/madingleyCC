#ifndef TMETABOLISMHETEROTROPH_H
#define TMETABOLISMHETEROTROPH_H
/** \file TMetabolismHeterotroph.h
 * \brief the TMetabolismHeterotroph header file
 */

//
//namespace Madingley
//{
/** \brief A formulation of the metabolism process 
\remarks Functional form and parameters taken from fitted relationship in Brown's (2004) Metabolic Theory of Ecology.
         Currently mass assigned to reproductive potential is not metabolised */
class MetabolismHeterotroph : public IMetabolismImplementation
    {
/** \brief Time units associated with the formulation of metabolism */
     StringProperty TimeUnitImplementation;
/** \brief Exponent describing the mass-dependency of metabolic rate */
        double MetabolismMassExponent;
/** \brief Normalization constatnt for metabolic rate  (independent of mass and temperature) */
        double NormalizationConstant;
/** \brief The activation energy of metabolism */
        double ActivationEnergy;
/** \brief Boltzmann's constant */
        double BoltzmannConstant;
/** \brief Scalar to convert energy in kJ to energy in grams mass */
        double EnergyScalar;

/** \brief Scalar to convert from the time step units used by this formulation of dispersal to global model time step units */
           DoubleProperty DeltaT;
/** \brief Include Utility class */
           UtilityFunctions Utilities;
/** \brief Constant to convert temperature in degrees Celsius to temperature in Kelvin */
          double TemperatureUnitsConvert;

/** \brief
Constructor for metabolism: assigns all parameter values
*/
MetabolismHeterotroph(string globalModelTimeStepUnit)
        {
           // Initialise ecological parameters for metabolism
           InitialiseMetabolismParameters();
           // Calculate the scalar to convert from the time step units used by this implementation of metabolism to the global  model time step units
           DeltaT = Utilities.ConvertTimeUnits(globalModelTimeStepUnit, TimeUnitImplementation());
           
           // Set the constant to convert temperature in degrees Celsius to Kelvin
           TemperatureUnitsConvert = 273.0;
        }

/** \brief Run metabolism for the acting cohort
@param gridCellCohorts The cohorts in the current grid cell 
@param gridCellStocks The stocks in the current grid cell 
@param actingCohort The position of the acting cohort in the jagged array of grid cell cohorts 
@param cellEnvironment The environment in the current grid cell 
@param deltas The sorted list to track changes in biomass and abundance of the acting cohort in this grid cell 
@param madingleyCohortDefinitions The definitions for cohort functional groups in the model 
@param madingleyStockDefinitions The definitions for the stock functional groups in the model 
@param currentTimestep The current model time step */
        void RunMetabolism(GridCellCohortHandler gridCellCohorts, GridCellStockHandler gridCellStocks, 
           vector<int> actingCohort, map<string, vector<double> > cellEnvironment, map<string, map<string, double>> 
           deltas, FunctionalGroupDefinitions madingleyCohortDefinitions, FunctionalGroupDefinitions madingleyStockDefinitions, 
           unsigned currentTimestep, unsigned currentMonth)
       {
           // Calculate metabolic loss for an individual and add the value to the delta biomass for metabolism
           deltas["biomass"]["metabolism"] = -CalculateIndividualMetabolicRate(gridCellCohorts[actingCohort].IndividualBodyMass(),
               cellEnvironment["Temperature"][currentMonth] + TemperatureUnitsConvert) * DeltaT();

           // If metabolic loss is greater than individual body mass after herbivory and predation, then set equal to individual body mass
           deltas["biomass"]["metabolism"] = max(deltas["biomass"]["metabolism"],-(gridCellCohorts[actingCohort].IndividualBodyMass() + deltas["biomass"]["predation"] + deltas["biomass"]["herbivory"]));

           // Add total metabolic loss for all individuals in the cohort to delta biomass for metabolism in the respiratory CO2 pool
           deltas["respiratoryCO2pool"]["metabolism"] = -deltas["biomass"]["metabolism"] * gridCellCohorts[actingCohort].CohortAbundance;
       }
/** \brief Initialises values for all ecological parameters for metabolism
\remarks Most parameters currently drawn from Brown's (2004) Metabolic Theory of Ecology
The scalar to convert kJ to grams mass currently a very rough estimate based on the calorific values
of fat, protein and carbohydrate*/
void InitialiseMetabolismParameters()
        {
           TimeUnitImplementation = "second";
           
           // Paramters from Brown's (2004) Metabolic Theory of Ecology
           MetabolismMassExponent = 0.71;
           NormalizationConstant = exp(20);
           ActivationEnergy = 0.69;
           BoltzmannConstant = 8.617e-5;

           // Currently a very rough estimate based on calorific values of fat, protein and carbohydrate
           EnergyScalar = 1.0 / 20000.0;

        }

/** \brief Calculate metabolic loss in grams for an individual
@param individualBodyMass The body mass of individuals in the acting cohort 
@param temperature The ambient temperature, in degrees Kelvin 
@return The metabolic loss for an individual*/
double CalculateIndividualMetabolicRate(double individualBodyMass, double temperature)
        {
           // Calculate metabolic loss in kJ
           double MetabolicLosskJ = NormalizationConstant * pow(individualBodyMass, MetabolismMassExponent) *
               exp(-(ActivationEnergy / (BoltzmannConstant * temperature)));

           // Return metabolic loss in grams
           return MetabolicLosskJ * EnergyScalar;

       }
   

       
    };
//}
#endif