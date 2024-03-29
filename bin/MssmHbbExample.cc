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
    mssmhbb.jetHistograms(config->nJetsMin(),"jet_matching");
    mssmhbb.jetHistograms(config->nJetsMin(),"bjet_matching3");
    mssmhbb.jetHistograms(config->nJetsMin(),"bjet_matching2");
    mssmhbb.add1DHistogram("jet_matching","matched_ranks","matched ranks",10,0,10);
    
    
    if ( wf == 1 )
    {
       for ( int i = 0 ; i < mssmhbb.nEvents() ; ++i )
        {
            if ( ! mssmhbb.event(i)               )  continue;    // read event, run selection/json
            if ( ! mssmhbb.triggerSelection()     )  continue;
            if ( ! mssmhbb.selectionJetId()       )  continue;
            if ( ! mssmhbb.selectionJetPileupId() )  continue;
            if ( ! mssmhbb.selectionNJets()       )  continue;
            if ( ! mssmhbb.jetCorrections()       )  continue;    // jet corrections
            if ( ! mssmhbb.selectionJet(1)        )  continue;    // jet selection
            if ( ! mssmhbb.selectionJet(2)        )  continue;    // jet selection
            if ( ! mssmhbb.selectionJetDr(1,2)    )  continue;    // jet deltaR selection
            if ( ! mssmhbb.selectionJetDeta(1,2)  )  continue;    // jet deltaEta selection
            if ( ! mssmhbb.onlineJetMatching(1)   )  continue;   // jet trg matching
            if ( ! mssmhbb.onlineJetMatching(2)   )  continue;   // jet trg matching
            mssmhbb.fillJetHistograms("kinematics");
            if ( ! mssmhbb.selectionBJet(1)       )  continue;   mssmhbb.actionApplyBtagSF(1);  // btag selection
            if ( ! mssmhbb.selectionBJet(2)       )  continue;   mssmhbb.actionApplyBtagSF(2);  // btag selection
            if ( ! mssmhbb.selectionJet(3)        )  continue;    // jet selection
            if ( ! mssmhbb.selectionJetDr(1,3)    )  continue;    // jet deltaR selection
            if ( ! mssmhbb.selectionJetDr(2,3)    )  continue;    // jet deltaR selection
            if ( ! mssmhbb.selectionBJet(3)       )  continue;   mssmhbb.actionApplyBtagSF(3); // btag selection
            mssmhbb.fillJetHistograms("btagging");
            if ( ! mssmhbb.onlineBJetMatching({1,2,3},2)  )  continue;   // bjet trg matching (more inclusive)
            if ( ! mssmhbb.selectionDiJetMass(1,2) ) continue;
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
            
            mssmhbb.fillJetHistograms("bjet_matching");
            mssmhbb.endSelection();
        } // end loop workflow 2
    } // end workflow 2
    if ( wf == 3 )
    {
       for ( int i = 0 ; i < mssmhbb.nEvents() ; ++i )
        {
            if ( ! mssmhbb.event(i)               )  continue;    // read event, run selection/json
            if ( ! mssmhbb.triggerSelection()     )  continue;
            if ( ! mssmhbb.selectionJetId()       )  continue;
            if ( ! mssmhbb.selectionJetPileupId() )  continue;
            if ( ! mssmhbb.selectionNJets()       )  continue;
            if ( ! mssmhbb.jetCorrections()       )  continue;    // jet corrections
            if ( ! mssmhbb.selectionJet(1)        )  continue;    // jet selection
            if ( ! mssmhbb.selectionJet(2)        )  continue;    // jet selection
            if ( ! mssmhbb.selectionJetDr(1,2)    )  continue;    // jet deltaR selection
            if ( ! mssmhbb.selectionJetDeta(1,2)  )  continue;    // jet deltaEta selection
            if ( ! mssmhbb.onlineJetMatching(1)   )  continue;   // jet trg matching
            if ( ! mssmhbb.onlineJetMatching(2)   )  continue;   // jet trg matching
            mssmhbb.fillJetHistograms("kinematics");
            if ( ! mssmhbb.selectionBJet(1)       )  continue;   mssmhbb.actionApplyBtagSF(1);  // btag selection
            if ( ! mssmhbb.selectionBJet(2)       )  continue;   mssmhbb.actionApplyBtagSF(2);  // btag selection
            mssmhbb.sortedBTagScorePlus3Jets();
            if ( ! mssmhbb.selectionJet(3)        )  continue;    // jet selection
            if ( ! mssmhbb.selectionJetDr(1,3)    )  continue;    // jet deltaR selection
            if ( ! mssmhbb.selectionJetDr(2,3)    )  continue;    // jet deltaR selection
            if ( ! mssmhbb.selectionBJet(3)       )  continue;   mssmhbb.actionApplyBtagSF(3); // btag selection
            mssmhbb.fillJetHistograms("btagging");
            if ( ! mssmhbb.onlineBJetMatching({1,2,3},2)  )  continue;   // bjet trg matching (more inclusive)
            if ( ! mssmhbb.selectionDiJetMass(1,2) ) continue;
            mssmhbb.endSelection();
        } // end loop workflow 3
    }
    if ( wf == 4 )
    {
       for ( int i = 0 ; i < mssmhbb.nEvents() ; ++i )
        {
            if ( ! mssmhbb.event(i)                       )  continue;    // read event, run selection/json
            if ( ! mssmhbb.triggerSelection()             )  continue;
            if ( ! mssmhbb.selectionPrimaryVertex()       )  continue;
            if ( ! mssmhbb.selectionJetId()               )  continue;
            if ( ! mssmhbb.selectionJetPileupId()         )  continue;
            if ( ! mssmhbb.selectionNJets()               )  continue;
            if ( ! mssmhbb.jetCorrections()               )  continue;    // jet corrections
            if ( ! mssmhbb.selectionJet(1)                )  continue;    // jet selection
            if ( ! mssmhbb.selectionJet(2)                )  continue;    // jet selection
            if ( ! mssmhbb.selectionJetDr(1,2)            )  continue;    // jet deltaR selection
            if ( ! mssmhbb.selectionJetDeta(1,2)          )  continue;    // jet deltaEta selection
            if ( ! mssmhbb.onlineJetMatching(1)           )  continue;   // jet trg matching
            if ( ! mssmhbb.onlineJetMatching(2)           )  continue;   // jet trg matching
            mssmhbb.fillJetHistograms("kinematics");      
            if ( ! mssmhbb.selectionBJet(1)               )  continue;   mssmhbb.actionApplyBtagSF(1);  // btag selection
            if ( ! mssmhbb.selectionBJet(2)               )  continue;   mssmhbb.actionApplyBtagSF(2);  // btag selection
            mssmhbb.sortedBTagScorePlus3Jets();
            mssmhbb.fsrCorrections(mssmhbb.mainJets(),mssmhbb.fsrCandidates());
            if ( ! mssmhbb.selectionJet(3)                )  continue;    // jet selection
            if ( ! mssmhbb.selectionJetDr(1,3)            )  continue;    // jet deltaR selection
            if ( ! mssmhbb.selectionJetDr(2,3)            )  continue;    // jet deltaR selection
            if ( ! mssmhbb.selectionBJet(3)               )  continue;   mssmhbb.actionApplyBtagSF(3); // btag selection
            mssmhbb.fillJetHistograms("btagging");
            if ( ! mssmhbb.onlineBJetMatching({1,2,3},2)  )  continue;   // bjet trg matching (more inclusive)
            if ( ! mssmhbb.selectionDiJetMass(1,2)        ) continue;
            mssmhbb.endSelection();
        } // end loop workflow 4
    }
    if ( wf == 5 )
    {
       for ( int i = 0 ; i < mssmhbb.nEvents() ; ++i )
        {
            if ( ! mssmhbb.event(i)                       )  continue;    // read event, run selection/json
            if ( ! mssmhbb.triggerSelection()             )  continue;
            if ( ! mssmhbb.selectionPrimaryVertex()       )  continue;
            if ( ! mssmhbb.selectionJetId()               )  continue;
            if ( ! mssmhbb.selectionJetPileupId()         )  continue;
            if ( ! mssmhbb.selectionNJets()               )  continue;
            if ( ! mssmhbb.jetCorrections()               )  continue;    // jet corrections
            if ( ! mssmhbb.selectionJet(1)                )  continue;    // jet selection
            if ( ! mssmhbb.selectionJet(2)                )  continue;    // jet selection
            if ( ! mssmhbb.selectionJetDr(1,2)            )  continue;    // jet deltaR selection
            if ( ! mssmhbb.selectionJetDeta(1,2)          )  continue;    // jet deltaEta selection
            if ( ! mssmhbb.onlineJetMatching(1)           )  continue;   // jet trg matching
            if ( ! mssmhbb.onlineJetMatching(2)           )  continue;   // jet trg matching
            mssmhbb.fillJetHistograms("kinematics");      
            if ( ! mssmhbb.selectionBJet(1)               )  continue;   mssmhbb.actionApplyBtagSF(1);  // btag selection
            if ( ! mssmhbb.selectionBJet(2)               )  continue;   mssmhbb.actionApplyBtagSF(2);  // btag selection
            mssmhbb.sortedBTagScorePlus3Jets();
            if ( ! mssmhbb.onlineBJetMatching({1,2,3},2)  )  continue;   // bjet trg matching (more inclusive)
            mssmhbb.fsrCorrections(mssmhbb.mainJets(),mssmhbb.fsrCandidates());
            if ( ! mssmhbb.selectionJet(3)                )  continue;    // jet selection
            if ( ! mssmhbb.selectionJetDr(1,3)            )  continue;    // jet deltaR selection
            if ( ! mssmhbb.selectionJetDr(2,3)            )  continue;    // jet deltaR selection
            if ( ! mssmhbb.selectionBJet(3)               )  continue;   mssmhbb.actionApplyBtagSF(3); // btag selection
            mssmhbb.fillJetHistograms("btagging");
            if ( ! mssmhbb.selectionDiJetMass(1,2)        ) continue;
            mssmhbb.endSelection();
        } // end loop workflow 5
    }
    if ( wf == 6 ) // btag weights
    {
        for ( int i = 0 ; i < mssmhbb.nEvents() ; ++i )
        {
            if ( ! mssmhbb.event(i)                       )  continue;    // read event, run selection/json
            if ( ! mssmhbb.triggerSelection()             )  continue;
            if ( ! mssmhbb.selectionPrimaryVertex()       )  continue;
            if ( ! mssmhbb.selectionJetId()               )  continue;
            if ( ! mssmhbb.selectionJetPileupId()         )  continue;
            if ( ! mssmhbb.selectionNJets()               )  continue;
            if ( ! mssmhbb.jetCorrections()               )  continue;    // jet corrections
            if ( ! mssmhbb.selectionJet(1)                )  continue;    // jet selection
            if ( ! mssmhbb.selectionJet(2)                )  continue;    // jet selection
            if ( ! mssmhbb.selectionJetDr(1,2)            )  continue;    // jet deltaR selection
            if ( ! mssmhbb.selectionJetDeta(1,2)          )  continue;    // jet deltaEta selection
            if ( ! mssmhbb.onlineJetMatching(1)           )  continue;   // jet trg matching
            if ( ! mssmhbb.onlineJetMatching(2)           )  continue;   // jet trg matching
            mssmhbb.fillJetHistograms("kinematics"); 
            if ( config->signalRegion() ) mssmhbb.sortedBTagScorePlus3Jets();
//            if ( fabs(mssmhbb.selectedJets()[2]->eta()) > 2.2 ) continue;
            if ( ! mssmhbb.selectionJet(3)                )  continue;    // jet selection
            if ( ! mssmhbb.onlineBJetMatching({1,2,3},2)  )  continue;   // bjet trg matching (more inclusive)
            mssmhbb.btagEfficiencyWeight();
            if ( ! config->signalRegion() ) // in CR apply reverse btag cut on the 3rd jet, otherwise already btag-weighted
                if ( ! mssmhbb.selectionBJet(3)               )  continue;
            mssmhbb.fsrCorrections(mssmhbb.mainJets(),mssmhbb.fsrCandidates());
            if ( ! mssmhbb.selectionJetDr(1,3)            )  continue;    // jet deltaR selection
            if ( ! mssmhbb.selectionJetDr(2,3)            )  continue;    // jet deltaR selection
            mssmhbb.fillJetHistograms("btagging");
            if ( ! mssmhbb.selectionDiJetMass(1,2)        ) continue;
            mssmhbb.endSelection();
        } // end loop workflow 6
    }
    if ( wf == 7 ) // btag cut w/ similar workflow to btag-weighted workflow 6
    {
        for ( int i = 0 ; i < mssmhbb.nEvents() ; ++i )
        {
            if ( ! mssmhbb.event(i)                       )  continue;    // read event, run selection/json
            if ( ! mssmhbb.triggerSelection()             )  continue;
            if ( ! mssmhbb.selectionPrimaryVertex()       )  continue;
            if ( ! mssmhbb.selectionJetId()               )  continue;
            if ( ! mssmhbb.selectionJetPileupId()         )  continue;
            if ( ! mssmhbb.selectionNJets()               )  continue;
            if ( ! mssmhbb.jetCorrections()               )  continue;    // jet corrections
            if ( ! mssmhbb.selectionJet(1)                )  continue;    // jet selection
            if ( ! mssmhbb.selectionJet(2)                )  continue;    // jet selection
            if ( ! mssmhbb.selectionJetDr(1,2)            )  continue;    // jet deltaR selection
            if ( ! mssmhbb.selectionJetDeta(1,2)          )  continue;    // jet deltaEta selection
            if ( ! mssmhbb.onlineJetMatching(1)           )  continue;   // jet trg matching
            if ( ! mssmhbb.onlineJetMatching(2)           )  continue;   // jet trg matching
            if ( config->signalRegion() ) mssmhbb.sortedBTagScorePlus3Jets();
            if ( ! mssmhbb.selectionJet(3)                )  continue;    // jet selection
//            if ( fabs(mssmhbb.selectedJets()[2]->eta()) > 2.2 ) continue;
            if ( ! mssmhbb.onlineBJetMatching({1,2,3},2)  )  continue;   // bjet trg matching (more inclusive)
            if ( ! mssmhbb.selectionBJet(1)               )  continue;   mssmhbb.actionApplyBtagSF(1);  // btag selection
            if ( ! mssmhbb.selectionBJet(2)               )  continue;   mssmhbb.actionApplyBtagSF(2);  // btag selection
            if ( ! mssmhbb.selectionBJet(3)               )  continue;   mssmhbb.actionApplyBtagSF(3);
            mssmhbb.fsrCorrections(mssmhbb.mainJets(),mssmhbb.fsrCandidates());
            if ( ! mssmhbb.selectionJetDr(1,3)            )  continue;    // jet deltaR selection
            if ( ! mssmhbb.selectionJetDr(2,3)            )  continue;    // jet deltaR selection
            mssmhbb.fillJetHistograms("btagging");
            if ( ! mssmhbb.selectionDiJetMass(1,2)        ) continue;
            mssmhbb.endSelection();
        } // end loop workflow 7
    }
    if ( wf == 8 ) // trigger selection 1,2 online btag - btag cut
    {
        for ( int i = 0 ; i < mssmhbb.nEvents() ; ++i )
        {
            if ( ! mssmhbb.event(i)                       )  continue;    // read event, run selection/json
            if ( ! mssmhbb.triggerSelection()             )  continue;
            if ( ! mssmhbb.selectionPrimaryVertex()       )  continue;
            if ( ! mssmhbb.selectionJetId()               )  continue;
            if ( ! mssmhbb.selectionJetPileupId()         )  continue;
            if ( ! mssmhbb.selectionNJets()               )  continue;
            if ( ! mssmhbb.jetCorrections()               )  continue;    // jet corrections
            if ( ! mssmhbb.selectionJet(1)                )  continue;    // jet selection
            if ( ! mssmhbb.selectionJet(2)                )  continue;    // jet selection
            if ( ! mssmhbb.selectionJetDr(1,2)            )  continue;    // jet deltaR selection
            if ( ! mssmhbb.selectionJetDeta(1,2)          )  continue;    // jet deltaEta selection
            if ( ! mssmhbb.onlineJetMatching(1)           )  continue;   // jet trg matching
            if ( ! mssmhbb.onlineJetMatching(2)           )  continue;   // jet trg matching
            if ( config->signalRegion() ) mssmhbb.sortedBTagScorePlus3Jets();
            if ( ! mssmhbb.selectionJet(3)                )  continue;    // jet selection
            if ( ! mssmhbb.onlineBJetMatching({1,2,3},3)  )  continue;             
            if ( ! mssmhbb.selectionBJet(1)               )  continue;   mssmhbb.actionApplyBtagSF(1);  // btag selection
            if ( ! mssmhbb.selectionBJet(2)               )  continue;   mssmhbb.actionApplyBtagSF(2);  // btag selection
            if ( ! mssmhbb.selectionBJet(3)               )  continue;   mssmhbb.actionApplyBtagSF(3);
            mssmhbb.fsrCorrections(mssmhbb.mainJets(),mssmhbb.fsrCandidates());
            if ( ! mssmhbb.selectionJetDr(1,3)            )  continue;    // jet deltaR selection
            if ( ! mssmhbb.selectionJetDr(2,3)            )  continue;    // jet deltaR selection
            if ( ! mssmhbb.selectionDiJetMass(1,2)        ) continue;
            mssmhbb.endSelection();
        } // end loop workflow 8
    }
    if ( wf == 9 ) // trigger selection 1,2 online btag - btag weighted
    {
        for ( int i = 0 ; i < mssmhbb.nEvents() ; ++i )
        {
            if ( ! mssmhbb.event(i)                       )  continue;    // read event, run selection/json
            if ( ! mssmhbb.triggerSelection()             )  continue;
            if ( ! mssmhbb.selectionPrimaryVertex()       )  continue;
            if ( ! mssmhbb.selectionJetId()               )  continue;
            if ( ! mssmhbb.selectionJetPileupId()         )  continue;
            if ( ! mssmhbb.selectionNJets()               )  continue;
            if ( ! mssmhbb.jetCorrections()               )  continue;    // jet corrections
            if ( ! mssmhbb.selectionJet(1)                )  continue;    // jet selection
            if ( ! mssmhbb.selectionJet(2)                )  continue;    // jet selection
            if ( ! mssmhbb.selectionJetDr(1,2)            )  continue;    // jet deltaR selection
            if ( ! mssmhbb.selectionJetDeta(1,2)          )  continue;    // jet deltaEta selection
            if ( ! mssmhbb.onlineJetMatching(1)           )  continue;   // jet trg matching
            if ( ! mssmhbb.onlineJetMatching(2)           )  continue;   // jet trg matching
            if ( config->signalRegion() ) mssmhbb.sortedBTagScorePlus3Jets();
            if ( ! mssmhbb.selectionJet(3)                )  continue;    // jet selection
            if ( ! mssmhbb.onlineBJetMatching({1,2,3},3)  )  continue;
            mssmhbb.actionApplyBtagEfficiency(1,1);
            mssmhbb.actionApplyBtagEfficiency(2,1);
            if ( config->signalRegion() ) 
            {
                mssmhbb.actionApplyBtagEfficiency(3,1);
            }
            else
            {
                if ( ! mssmhbb.selectionBJet(3)               )  continue;   mssmhbb.actionApplyBtagSF(3);
            }
            mssmhbb.fsrCorrections(mssmhbb.mainJets(),mssmhbb.fsrCandidates());
            if ( ! mssmhbb.selectionJetDr(1,3)            )  continue;    // jet deltaR selection
            if ( ! mssmhbb.selectionJetDr(2,3)            )  continue;    // jet deltaR selection
            if ( ! mssmhbb.selectionDiJetMass(1,2)        ) continue;
            mssmhbb.endSelection();
        } // end loop workflow 9
    }

} 

// end main
     



