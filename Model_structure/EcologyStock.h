#ifndef ECOLOGYSTOCK_H
#define ECOLOGYSTOCK_H
/** \file EcologyStock.h
 * \brief the EcologyStock header file
 */






//

//
//namespace Madingley
//{
/** \brief
//    /// A class to specify, initialise and run ecological processes pertaining to stocks
//    /// </summary>
//    class EcologyStock
//    {
/** \brief
An instance of the Autotroph Processor for this model
*/
//        AutotrophProcessor MarineNPPtoAutotrophStock;
//
/** \brief
An instance of the plant model class
*/
//        RevisedTerrestrialPlantModel DynamicPlantModel;
//
/** \brief
An instance of the class for human appropriation of NPP
*/
//        HumanAutotrophMatterAppropriation HANPP;
//
//
//         void InitializeEcology()
//        {
//            //Initialize the autotrophprocessor
//            MarineNPPtoAutotrophStock = new AutotrophProcessor();
//
//            // Initialise the plant model
//            DynamicPlantModel = new RevisedTerrestrialPlantModel();
//
//            // Initialise the human NPP appropriation class
//            HANPP = new HumanAutotrophMatterAppropriation();
//
//        }
//
//
/** \brief
Run ecological processes that operate on stocks within a single grid cell
*/
//        ///@param gridCellStocks The stocks in the current grid cell 
//        ///@param actingStock The acting stock 
//        ///@param cellEnvironment The stocks in the current grid cell 
//        ///@param environmentalDataUnits List of units associated with the environmental variables 
//        ///@param humanNPPExtraction Name of the human appropriation of NPP scenario to run 
//        ///@param madingleyStockDefinitions The definitions for stock functional groups in the model 
//        ///@param currentTimeStep The current model time step 
//        ///@param globalModelTimeStepUnit The time step unit used in the model 
//        ///@param trackProcesses Whether to track properties of ecological processes 
//        ///@param tracker An instance of the ecological process tracker 
//        ///@param currentMonth The current model month 
//         void RunWithinCellEcology(GridCellStockHandler gridCellStocks, int[] actingStock, map<string, double[]> cellEnvironment, 
//            map<string,string> environmentalDataUnits, string humanNPPExtraction, FunctionalGroupDefinitions madingleyStockDefinitions, 
//            unsigned currentTimeStep, string globalModelTimeStepUnit, Boolean trackProcesses, ProcessTracker tracker, GlobalProcessTracker globalTracker, unsigned currentMonth, 
//            string outputDetail, bool specificLocations)
//        {
//            if (madingleyStockDefinitions.GetTraitNames("Realm", actingStock[0]) == "marine")
//            {
//                // Run the autotroph processor
//                MarineNPPtoAutotrophStock.ConvertNPPToAutotroph(cellEnvironment, gridCellStocks, actingStock, environmentalDataUnits["LandNPP"], 
//                    environmentalDataUnits["OceanNPP"], currentTimeStep,globalModelTimeStepUnit,tracker,globalTracker ,outputDetail,specificLocations,currentMonth);
//            }
//            else if (madingleyStockDefinitions.GetTraitNames("Realm", actingStock[0]) == "terrestrial")
//            {
//                // Run the dynamic plant model to update the leaf stock for this time step
//                DynamicPlantModel.UpdateLeafStock(cellEnvironment, gridCellStocks, actingStock, currentTimeStep, madingleyStockDefinitions.
//                    GetTraitNames("leaf strategy", actingStock[0]).Equals("deciduous"), globalModelTimeStepUnit, tracker, globalTracker, currentMonth,
//                    outputDetail, specificLocations);
//
//                // Apply human appropriation of NPP
//                HANPP.RemoveHumanAppropriatedMatter(cellEnvironment, humanNPPExtraction, gridCellStocks, actingStock, currentTimeStep, currentMonth);
//
//            }
//            else
//            {
//                Debug.Fail("Stock must be classified as belonging to either the marine or terrestrial realm");
//            }
//        }
//    }
//}
#endif
