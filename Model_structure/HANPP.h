#ifndef HANPP_H
#define HANPP_H
/** \file HANPP.h
 * \brief the HANPP header file
 */






//
//namespace Madingley
//{
/** \brief
//    /// Removes autotroph matter appropriated by humans from a grid cell's autotroph stocks
//    /// </summary>
//    /// <remarks>Assumes that autotroph matter is appropriated evenly from different stocks in proportion to their biomass</remarks>
//     class HumanAutotrophMatterAppropriation
//    {
/** \brief
Constructor for human appropriation of autotroph matter
*/
//         HumanAutotrophMatterAppropriation()
//        {
//        }
//
/** \brief
Remove human appropriated matter from the grid cell autotroph stocks
*/
@param cellEnvironment The environment in the current grid cell 
@param humanNPPExtraction The type of NPP extraction to apply: 'no' = no removal; 'hanpp' = appropriated NPP estimate from input map; or proportion of total NPP 
@param gridCellStocks The stocks in the current grid cell 
@param actingStock The position of the acting stock in the jagged array of grid cell stocks 
@param currentTimestep The current model time step 
//         void RemoveHumanAppropriatedMatter(map<string,double[]> cellEnvironment, string humanNPPExtraction, GridCellStockHandler 
//            gridCellStocks, int[] actingStock, unsigned currentTimestep,unsigned currentMonth)
//        {
//            // Factor to convert NPP from units per m2 to units per km2
//            double m2Tokm2Conversion = 1000000.0;
//
//
//            if (humanNPPExtraction == "hanpp")
//            {
//                // Loop over stocks in the grid cell and calculate the total biomass of all stocks
//                double TotalAutotrophBiomass = 0.0;
//                foreach (var stockFunctionalGroup in gridCellStocks)
//                {
//                    for (int i = 0; i < stockFunctionalGroup.Count; i++)
//                    {
//                        TotalAutotrophBiomass += stockFunctionalGroup[i].TotalBiomass;
//                    }
//                }
//
//                // Get the total amount of NPP appropriated by humans from this cell
//                double HANPP = cellEnvironment["HANPP"][0] * cellEnvironment["Seasonality"][currentMonth];
//
//                // If HANPP value is missing, then assume zero
//                if (HANPP == cellEnvironment["Missing Value"][0]) HANPP = 0.0;
//
//                // Allocate HANPP for this stock according to the proportion of total autotroph biomass that the stock represents
//                if (TotalAutotrophBiomass == 0.0)
//                {
//                    HANPP = 0.0;
//                }
//                else
//                {
//                    HANPP *= (gridCellStocks[actingStock].TotalBiomass / TotalAutotrophBiomass);
//                }
//
//
//                // Convert gC/m2/month to gC/km2/month
//                HANPP *= m2Tokm2Conversion;
//
//                // Multiply by cell area (in km2) to get g/cell/day
//                HANPP *= cellEnvironment["Cell Area"][0];
//
//                // Convert from gC to g dry matter
//                double DryMatterAppropriated = HANPP * 2;
//
//                // Convert from g dry matter to g wet matter
//                double WetMatterAppropriated = DryMatterAppropriated * 2;
//
//
//                // Remove human appropriated NPP from total NPP, reduce NPP to 10% available to herbivores and then add to autotroph biomass
//                gridCellStocks[actingStock].TotalBiomass -= WetMatterAppropriated;
//
//                if (gridCellStocks[actingStock].TotalBiomass < 0.0) gridCellStocks[actingStock].TotalBiomass = 0.0;
//
//            }
//            else if (humanNPPExtraction == "no")
//            {
//            }
//            else
//            {
//                Debug.Assert((Convert.ToDouble(humanNPPExtraction) >= 0.0) && (Convert.ToDouble(humanNPPExtraction) <= 1.0),
//                        "Human NPP extraction must be specified in the scenarios file as 'no', 'hanpp' or a proportion between zero and one");
//
//                // Get the proportion of plant matter appropriated as specified
//                double ProportionMatterAppropriated = Convert.ToDouble(humanNPPExtraction);
//
//                // Get the absolute amount of NPP appropriated based on this
//                double MatterAppropriated = gridCellStocks[actingStock].TotalBiomass * ProportionMatterAppropriated;
//
//                // Remove human appropriated NPP from total NPP, reduce NPP to 10% available to herbivores and then add to autotroph biomass
//                gridCellStocks[actingStock].TotalBiomass -= MatterAppropriated;
//            }
//            
//
//        }
//        
//
//    }
//}
#endif
