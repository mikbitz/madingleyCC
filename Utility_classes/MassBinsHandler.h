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