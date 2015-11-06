/** \file Stock.h
 * \brief the Stock implementation file
 */
#include <RevisedTerrestrialPlantModel.h>
using namespace std;
    //----------------------------------------------------------------------------------------------
    //Methods
    //----------------------------------------------------------------------------------------------
    //Constructor
    Stock::Stock(FunctionalGroupDefinitions& StockDefinitions, const unsigned FunctionalGroup, map<string, vector<double>>&Environment, bool& success) {
        
            FunctionalGroupIndex = FunctionalGroup;

            // Get the individual body masses for organisms in each stock functional group
            IndividualBodyMass = StockDefinitions.GetBiologicalPropertyOneFunctionalGroup("individual mass",FunctionalGroup);

            double Mass = 0;
            success = false;
            // If it is a functional group that corresponds to the current realm, then seed the stock
            if (Environment["Realm"][0] == 1.0 && Environment["Precipitation"][0] != Environment["Missing Value"][0] && Environment["Temperature"][0] != Environment["Missing Value"][0]) {
                if (StockDefinitions.GetTraitNames("Realm", FunctionalGroup) == "terrestrial") {
                    // An instance of the terrestrial carbon model class
                    RevisedTerrestrialPlantModel PlantModel;

                    // Calculate predicted leaf mass at equilibrium for this stock
                    TotalBiomass = PlantModel.CalculateEquilibriumLeafMass(Environment, StockDefinitions.GetTraitNames("leaf strategy", FunctionalGroup) == "deciduous");
                    success = true;
                }
            } else if (Environment["Realm"][0] == 2.0 && Environment["NPP"][0] != Environment["Missing Value"][0]) {
                if (StockDefinitions.GetTraitNames("Realm", FunctionalGroup) == "marine") {
                    TotalBiomass = 1.e12;
                    success = true;
                }
            }
        }



