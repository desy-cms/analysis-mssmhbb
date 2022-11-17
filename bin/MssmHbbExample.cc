#include "Analysis/MssmHbb/interface/MssmHbbAnalyser.h"

using namespace analysis;
using namespace analysis::tools;
using namespace analysis::mssmhbb;

int main(int argc, char ** argv)
{
    TH1::SetDefaultSumw2();  // proper treatment of errors when scaling histograms

    MssmHbbAnalyser mssmhbb(argc,argv);
    
    auto config = mssmhbb.config();
    auto wf = config->workflow();

    mssmhbb.jetHistograms(config->nJetsMin(),"kinematics");
    mssmhbb.jetHistograms(config->nJetsMin(),"btagging");
    mssmhbb.jetHistograms(config->nJetsMin(),"trigger_matching");
    
    if ( wf == 1 )
    {
       for ( int i = 0 ; i < mssmhbb.nEvents() ; ++i )
        {
            if ( ! mssmhbb.event(i)               )  continue;    // read event, run selection/json
            if ( ! mssmhbb.triggerSelection()     )  continue;
            if ( ! mssmhbb.selectionJetId()       )  continue;
            if ( ! mssmhbb.selectionJetPileupId() )  continue;
            if ( ! mssmhbb.jetCorrections()       )  continue;    // jet corrections
            if ( ! mssmhbb.selectionJet(1)        )  continue;    // jet selection
            if ( ! mssmhbb.selectionJet(2)        )  continue;    // jet selection
            if ( ! mssmhbb.selectionJet(3)        )  continue;    // jet selection
            if ( ! mssmhbb.selectionJetDr(1,2)    )  continue;    // jet deltaR selection
            if ( ! mssmhbb.selectionJetDr(1,3)    )  continue;    // jet deltaR selection
            if ( ! mssmhbb.selectionJetDr(2,3)    )  continue;    // jet deltaR selection
            if ( ! mssmhbb.selectionJetDeta(1,2)  )  continue;    // jet deltaEta selection
            mssmhbb.fillJetHistograms("kinematics");
            if ( ! mssmhbb.selectionBJet(1)       )  continue;   mssmhbb.actionApplyBtagSF(1);  // btag selection
            if ( ! mssmhbb.selectionBJet(2)       )  continue;   mssmhbb.actionApplyBtagSF(2);  // btag selection
            if ( ! mssmhbb.selectionBJet(3)       )  continue;   mssmhbb.actionApplyBtagSF(3); // btag selection
            mssmhbb.fillJetHistograms("btagging");
            if ( ! mssmhbb.onlineJetMatching(1)   )  continue;   // jet trg matching
            if ( ! mssmhbb.onlineJetMatching(2)   )  continue;   // jet trg matching
            if ( ! mssmhbb.onlineBJetMatching({1,2,3},2)  )  continue;   // bjet trg matching
            
            mssmhbb.fillJetHistograms("trigger_matching");
            mssmhbb.endSelection();
        } // end loop workflow 1
    }

    if ( wf == 2 )
    {
       for ( int i = 0 ; i < mssmhbb.nEvents() ; ++i )
        {
            if ( ! mssmhbb.event(i)               )  continue;    // read event, run selection/json
            if ( ! mssmhbb.triggerSelection()     )  continue;
            if ( ! mssmhbb.selectionJetId()       )  continue;
            if ( ! mssmhbb.selectionJetPileupId() )  continue;
            if ( ! mssmhbb.jetCorrections()       )  continue;    // jet corrections
            if ( ! mssmhbb.selectionJet(1)        )  continue;    // jet selection
            if ( ! mssmhbb.selectionJet(2)        )  continue;    // jet selection
            if ( ! mssmhbb.selectionJet(3)        )  continue;    // jet selection
            if ( ! mssmhbb.selectionJetDr(1,2)    )  continue;    // jet deltaR selection
            if ( ! mssmhbb.selectionJetDr(1,3)    )  continue;    // jet deltaR selection
            if ( ! mssmhbb.selectionJetDr(2,3)    )  continue;    // jet deltaR selection
            if ( ! mssmhbb.selectionJetDeta(1,2)  )  continue;    // jet deltaEta selection
            mssmhbb.fillJetHistograms("kinematics");
            if ( ! mssmhbb.selectionBJet(1)       )  continue;   mssmhbb.actionApplyBtagSF(1);  // btag selection
            if ( ! mssmhbb.selectionBJet(2)       )  continue;   mssmhbb.actionApplyBtagSF(2);  // btag selection
            if ( ! mssmhbb.selectionBJet(3)       )  continue;   mssmhbb.actionApplyBtagSF(3); // btag selection
            mssmhbb.fillJetHistograms("btagging");
            if ( ! mssmhbb.onlineJetMatching(1)   )  continue;   // jet trg matching
            if ( ! mssmhbb.onlineJetMatching(2)   )  continue;   // jet trg matching
            if ( ! mssmhbb.onlineBJetMatching(1)  )  continue;   // bjet trg matching
            if ( ! mssmhbb.onlineBJetMatching(2)  )  continue;   // bjet trg matching
            
            mssmhbb.fillJetHistograms("trigger_matching");
            mssmhbb.endSelection();
        } // end loop workflow 2
    } // end workflow 2

} 

// end main
     



