#ifndef AUTOTROPHPROCESSOR_H
#define AUTOTROPHPROCESSOR_H
/** \file AutotrophProcessor.h
 * \brief the AutotrophProcessor header file
 */






//
//namespace Madingley
//{
/** \brief
//    /// Class for converting primary productivity estimates to autotroph biomass
//    /// </summary>
//     class AutotrophProcessor
//    {
/** \brief
Instance of the class to perform general functions
*/
//        private UtilityFunctions Utilities;
//
//        // Conversion ratio for phytoplankton from grams carbon to grams wet weight
/** \brief
Factor to convert phytoplankton biomass from grams carbon to grams wet weight
*/
<remarks>Currently derived from Ho et al. (2003) J. Phycol., Dalsgaard and Pauly (1997) and Strickland (1966)</remarks>
//        private const double _PhytoplanktonConversionRatio = 10;
/** \brief
Get the conversion ratio for phytoplankton from grams carbon to grams wet weight
*/
//         double PhytoplanktonConversionRatio { get {return _PhytoplanktonConversionRatio; } }
//
/** \brief
Factor to convert NPP from units per m^2 to units per km^2
*/
//        private const double _MsqToKmSqConversion = 1000000.0;
/** \brief
Get the factor to convert NPP from units per m^2 to units per km^2
*/
//         double MsqToKmSqConversion { get { return _MsqToKmSqConversion; } }
//
/** \brief
Constructor for the autotroph processor: initialises necessary classes
*/
//         AutotrophProcessor()
//        {
//            // Initialise the utility functions
//            Utilities = new UtilityFunctions();
//        }
//
/** \brief
Convert NPP estimate into biomass of an autotroph stock
*/
@param cellEnvironment The environment of the current grid cell 
@param gridCellStockHandler The stock handler for the current stock 
@param actingStock The location of the stock to add biomass to 
@param terrestrialNPPUnits The units of the terrestrial NPP data 
@param oceanicNPPUnits The units of the oceanic NPP data 
@param currentTimestep The current model time step 
@param GlobalModelTimeStepUnit The time step unit used in the model 
//         void ConvertNPPToAutotroph(map<string,double[]> cellEnvironment, GridCellStockHandler gridCellStockHandler, int[] 
//            actingStock, string terrestrialNPPUnits, string oceanicNPPUnits, unsigned currentTimestep, string GlobalModelTimeStepUnit,
//            ProcessTracker trackProcesses, GlobalProcessTracker globalTracker, string outputDetail, bool specificLocations,unsigned currentMonth)
//        {
//
//            // Get NPP from the cell environment
//            double NPP = cellEnvironment["NPP"][currentMonth];
//
//            // If NPP is a mssing value then set to zero
//            if (NPP == cellEnvironment["Missing Value"][0]) NPP = 0.0;
//
//            // Check that this is an ocean cell
//            if (cellEnvironment["Realm"][0] == 2.0)
//            {
//                // Check that the units of oceanic NPP are gC per m2 per day
//                Debug.Assert(oceanicNPPUnits == "gC/m2/day", "Oceanic NPP data are not in the correct units for this formulation of the model");
//
//                // Convert to g/cell/month
//                NPP *= _MsqToKmSqConversion;
//
//                // Multiply by cell area to get g/cell/day
//                NPP *= cellEnvironment["Cell Area"][0];
//
//                // Convert to g wet matter, assuming carbon content of phytoplankton is 10% of wet matter
//                NPP *= _PhytoplanktonConversionRatio;
//
//                // Finally convert to g/cell/month and add to the stock totalbiomass
//                NPP *= Utilities.ConvertTimeUnits(GlobalModelTimeStepUnit, "day");
//                gridCellStockHandler[actingStock].TotalBiomass += NPP;
//
//                if (trackProcesses.TrackProcesses && (outputDetail == "high") && specificLocations)
//                {
//                    trackProcesses.TrackPrimaryProductionTrophicFlow((unsigned)cellEnvironment["LatIndex"][0], (unsigned)cellEnvironment["LonIndex"][0],
//                        NPP);
//                }
//
//                if (globalTracker.TrackProcesses)
//                {
//                    globalTracker.RecordNPP((unsigned)cellEnvironment["LatIndex"][0], (unsigned)cellEnvironment["LonIndex"][0],
//                            NPP / cellEnvironment["Cell Area"][0]);
//                }
//
//                // If the biomass of the autotroph stock has been made less than zero (i.e. because of negative NPP) then reset to zero
//                if (gridCellStockHandler[actingStock].TotalBiomass < 0.0)
//                    gridCellStockHandler[actingStock].TotalBiomass = 0.0;
//            }
//            // Else if neither on land or in the ocean
//            else
//            {
//                Debug.Fail("This is not a marine cell!");
//                // Set the autotroph biomass to zero
//                gridCellStockHandler[actingStock].TotalBiomass = 0.0;
//            }
//            Debug.Assert(gridCellStockHandler[actingStock].TotalBiomass >= 0.0, "stock negative");
//        }
//
//    }
//}
#endif