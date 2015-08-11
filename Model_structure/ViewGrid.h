#ifndef VIEWGRID_H
#define VIEWGRID_H
/** \file ViewGrid.h
 * \brief the ViewGrid header file
 */





//




//
//namespace Madingley
//{
/** \brief
//    /// This class is for viewing gridded data such as environmental data layers or state variables; it pauses program execution while the viewer is open
//    /// </summary>
//     class ViewGrid
//    {
/** \brief
An instance of the class to convert data between arrays and SDS objects
*/
//        private ArraySDSConvert DataConverter;
//
/** \brief
Constructor for the grid viewer: initialses relevant objects
*/
//         ViewGrid()
//        {
//            DataConverter = new ArraySDSConvert();
//        }
//
/** \brief
Copy an georeferenced array (should be by reference!) to a grid in order to view it, then spawn the data set viewer
@param gridToView The grid to be viewed 
@param variableName The name of the variable to be viewed 
@param lats A vector of latitudes associated with the grid 
@param lons A vector of longitudes associated with the grid 
@param gridMissingValue The missing value for the grid to view 
*/
//         void PauseProgram(ref double[,] gridToView, string variableName, float[] lats, float[] lons, double gridMissingValue)
//        {
//            // Create a new data set, set it to commit changes manually
//            var DataSetToView = DataSet.Open("msds:memory");
//            DataSetToView.IsAutocommitEnabled = false;
//
//            // Convert the grid to an SDS structure
//            //ConvertSDSArray.GridToSDS(ref gridToView, variableName, lats, lons, gridMissingValue, ref DataSetToView, false);
//            DataConverter.Array2DToSDS2D(gridToView, variableName, lats, lons, gridMissingValue, DataSetToView);
//
//            // Open the viewer
//            DataSetToView.View();
//
//            // Remove the temporary data set from memory
//            DataSetToView.Dispose();
//        }
//
/** \brief
Provides a snapshot view of an SDS
*/
@param DataSetToView The name of the SDS to view 
@param handle An object handle for the viewer instance; send the same handle to prevent multiple instances of SDS viewer opening 
<todoD>Need to update to be able to select which variable to view</todoD>
<todoD>Pass sleep length</todoD>
//         void SnapshotView(ref DataSet DataSetToView, ref object handle)
//        {
//            // Open the snapshot viewer
//            handle = DataSetToView.ViewSnapshot("", handle);
//            
//            // Slow down computation
//            System.Threading.Thread.Sleep(250);
//
//        }
//
/** \brief
Asynchronously views an SDS
*/
@param DataSetToView The name of the SDS to view 
@param viewingParameters A string of viewing parameters ('hints') to pass to SDS viewer 
<todoD>Need to update to be able to select which variable to view</todoD>
<todoD>Pass sleep length</todoD>
<todoD>UPdate title on each timestep</todoD>
//         void AsynchronousView(ref DataSet DataSetToView, string viewingParameters)
//        {
//            DataSetToView.SpawnViewer(viewingParameters);
//            // Slow down computation
//            //System.Threading.Thread.Sleep(10);
//
//        }
//
/** \brief
Asynchronously views an SDS
*/
@param DataSetToView The name of the SDS to view 
<todoD>Need to update to be able to select which variable to view</todoD>
<todoD>Pass sleep length</todoD>
<todoD>UPdate title on each timestep</todoD>
//         void AsynchronousView(ref DataSet DataSetToView)
//        {
//            DataSetToView.SpawnViewer();
//            // Slow down computation
//            //System.Threading.Thread.Sleep(10);
//
//        }
//
//
//    }
//}
#endif
