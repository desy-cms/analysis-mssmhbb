/**\class MssmHbb MssmHbbAnalyser.cc Analysis/Tools/src/MssmHbbAnalyser.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Roberval Walsh Bastos Rangel
//         Created:  Mon, 20 Oct 2014 14:24:08 GMT
//
//

// system include files
#include <iostream>
// 
// user include files
#include "Analysis/MssmHbb/interface/MssmHbbAnalyser.h"


//
// class declaration
//

using namespace analysis;
using namespace analysis::tools;
using namespace analysis::mssmhbb;

//
// constructors and destructor
//
MssmHbbAnalyser::MssmHbbAnalyser()
{
}

MssmHbbAnalyser::MssmHbbAnalyser(int argc, char * argv[]) : Analyser(argc,argv)
{
   histograms("jet",njetsmin_);
   
}

MssmHbbAnalyser::~MssmHbbAnalyser()
{
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)
}


//
// member functions
//
// ------------ method called for each event  ------------


bool MssmHbbAnalyser::event(const int & i)
{
   if ( ! Analyser::event(i) ) return false;
   
   return true;
}


void MssmHbbAnalyser::histograms(const std::string & obj, const int & n)
{
   Analyser::histograms(obj,n);
   if ( obj == "jet" )
   {
   }
   
   
}
