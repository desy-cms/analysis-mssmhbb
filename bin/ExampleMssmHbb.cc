#include "Analysis/MssmHbb/interface/MssmHbbAnalyser.h"

// this example is a dijet analysis in which the leading jet
// is required to have a muon, and both of the leading
// jets are btagged

using namespace std;
using namespace analysis;
using namespace analysis::mssmhbb;

int main(int argc, char ** argv)
{
   TH1::SetDefaultSumw2();  // proper treatment of errors when scaling histograms
   
   MssmHbbAnalyser analyser(argc,argv);
   
// HISTOGRAMS   
   // create some predefined jet histograms
   analyser.jetHistograms(3,"initial");
   analyser.jetHistograms(3,"final");
   // create some predefined muon histograms
   // muon histograms still not available
   
   for ( int i = 0 ; i < analyser.nEvents() ; ++i )
   {
      if ( ! analyser.event(i)                  )   continue;
      
// JETS pre-selection 
      // jet identification selection
      if ( ! analyser.selectionJetId()          )   continue;
      if ( ! analyser.selectionJetPileupId()    )   continue;
      if ( ! analyser.selectionNJets()          )   continue;
      
      analyser.fillJetHistograms("initial");
     
// CORRECTIONS
   // b energy regression
      if ( analyser.config()->bRegression() )  analyser.actionApplyBjetRegression();
   // Jet energy resolution smearing
      analyser.actionApplyJER();
      
// MAIN SELECTION
   // JETS
      //  1st and 2nd jet kinematic selection, pt and eta
      if ( ! analyser.selectionJet(1)          )   continue;
      if ( ! analyser.selectionJet(2)          )   continue;
      if ( ! analyser.selectionJet(3)          )   continue;
      // delta R jet selection
      if ( ! analyser.selectionJetDr(1,2)      )   continue;
      if ( ! analyser.selectionJetDr(1,3)      )   continue;
      if ( ! analyser.selectionJetDr(2,3)      )   continue;
      if ( ! analyser.selectionJetDeta(1,2)    )   continue;
      // btag
      if ( ! analyser.selectionBJet(1)         )   continue;
      if ( ! analyser.selectionBJet(2)         )   continue;
      if ( ! analyser.selectionBJet(3)         )   continue;
// TRIGGER selection
      if ( ! analyser.selectionHLT()           )   continue;
      if ( ! analyser.selectionL1 ()           )   continue;  // to be used mainly in case of "OR" of seeds
      // jets 1 and 2 matching to online jets and btag
      if ( ! analyser.onlineJetMatching(1)     )   continue;
      if ( ! analyser.onlineJetMatching(2)     )   continue;
      if ( ! analyser.onlineBJetMatching(1)    )   continue;
      if ( ! analyser.onlineBJetMatching(2)    )   continue;
            
      // btag scale factors
      analyser.actionApplyBtagSF(1);
      analyser.actionApplyBtagSF(2);
      analyser.actionApplyBtagSF(3);

// HISTOGRAMS
      analyser.fillJetHistograms("final");
      
   }  //end event loop
} // end main
      
