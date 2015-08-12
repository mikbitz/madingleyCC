#ifndef MASSBINHANDLER_H
#define MASSBINHANDLER_H
#include <fstream>
//namespace Madingley
//{
/** \brief Handles the mass bins to be used in model outputs */

class MassBinsHandler
    {
    public:
/** \brief The number of mass bins to be used for outputs */
        int NumMassBins = 50;       
/** \brief A vector containing the masses correpsonding to the mass bins */
        vector<float> MassBins;
//
/** \brief Sets up mass bins based on an input file
@param massBinsFile The filename containing the mass bin information 
@param outputPath The path to the output folder to copy the mass bins definition file to
*/
       void SetUpMassBins(string massBinsFile, string outputPath)
       {
//            // Construct file name
//            string FileString = "msds:csv?file=input/Model setup/" + massBinsFile + "&openMode=readOnly";
//
//            //Copy the Mass bin definitions file to the output directory
//            if(System.IO.File.Exists(outputPath+ massBinsFile))
//                System.IO.File.Copy("input/Model setup/" + massBinsFile, outputPath + massBinsFile, true);
//
//            // Read in the data
//            DataSet InternalData = DataSet.Open(FileString);
//
//            //Copy the values for this variable into an array
//            var TempValues = InternalData.Variables[0].GetData();
//            NumMassBins = TempValues.Length;
//            MassBins = new float[TempValues.Length];
//
//            for (int i = 0; i < TempValues.Length; i++)
//            {
//                MassBins[i] = Convert.ToSingle(TempValues.GetValue(i));
//            }
           ifstream massFile(massBinsFile.c_str());
           string title;
           getline (massFile,title);
           float f;
           while (!massFile.eof()){
               massFile>>f;
               if (!massFile.eof())MassBins.push_back(f);
           }
           // Sort the array of mass bins
           sort(MassBins.begin(),MassBins.end());
       }

/** \brief Returns the mass bins copied from file
@return the mass bins copied from file
*/
       vector<float> const& GetSpecifiedMassBins() const
       {
           return MassBins;
       }

};
#endif