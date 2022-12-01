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
            mssmhbb.fillJetHistograms("jet_matching");
            if ( ! mssmhbb.selectionDiJetMass(1,2) ) continue;
            auto bmatched = mssmhbb.onlineBJetMatching({1,2,3});
            if ( bmatched.size() > 1 )
            {
                if ( bmatched[0] == 1 && bmatched[1] == 2 ) mssmhbb.fill1DHistogram("jet_matching","matched_ranks",0.,mssmhbb.weight());
                if ( bmatched[0] == 1 && bmatched[1] == 3 ) mssmhbb.fill1DHistogram("jet_matching","matched_ranks",1.,mssmhbb.weight());
                if ( bmatched[0] == 2 && bmatched[1] == 3 ) mssmhbb.fill1DHistogram("jet_matching","matched_ranks",2.,mssmhbb.weight());
            }
            if ( ! mssmhbb.onlineBJetMatching({1,2,3},2)  )  continue;   // bjet trg matching (more inclusive)
            mssmhbb.fillJetHistograms("bjet_matching3");
            if ( ! mssmhbb.onlineBJetMatching({1,2},2)  )  continue;   // bjet trg matching (exclusive)
            mssmhbb.fillJetHistograms("bjet_matching2");
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
            // auto seljets = mssmhbb.selectedJets();
            // for ( auto & jet : seljets )
            // {
            //     std::cout << "jet pt = " << jet->pt() << ", btag = " << jet->btag() << ", flavour = " << jet->flavour() << std::endl;
            // }
            // std::cout << "---------------" << std::endl;
            // auto two_leading_jets = mssmhbb.keepSelectedJets({1,2});
            // auto plus3_leading_jets = mssmhbb.removeSelectedJets({1,2});
            // auto btag_sorted_jets = mssmhbb.btagSortedJets(plus3_leading_jets);
            // auto new_selectedJets = mssmhbb.concatenateJets(two_leading_jets,btag_sorted_jets);
            mssmhbb.sortedBTagScorePlus3Jets();
            // auto newseljets = mssmhbb.selectedJets();
            // for ( auto & jet : newseljets )
            // {
            //     std::cout << "jet pt = " << jet->pt() << ", btag = " << jet->btag() << ", flavour = " << jet->flavour() << std::endl;
            // }
            // std::cout << "================== index = " << i << std::endl;
            if ( ! mssmhbb.selectionBJet(3)       )  continue;   mssmhbb.actionApplyBtagSF(3); // btag selection
            mssmhbb.fillJetHistograms("btagging");
            if ( ! mssmhbb.onlineJetMatching(1)   )  continue;   // jet trg matching
            if ( ! mssmhbb.onlineJetMatching(2)   )  continue;   // jet trg matching
            mssmhbb.fillJetHistograms("jet_matching");
            if ( ! mssmhbb.selectionDiJetMass(1,2) ) continue;
            auto bmatched = mssmhbb.onlineBJetMatching({1,2,3});
            if ( bmatched.size() > 1 )
            {
                if ( bmatched[0] == 1 && bmatched[1] == 2 ) mssmhbb.fill1DHistogram("jet_matching","matched_ranks",0.,mssmhbb.weight());
                if ( bmatched[0] == 1 && bmatched[1] == 3 ) mssmhbb.fill1DHistogram("jet_matching","matched_ranks",1.,mssmhbb.weight());
                if ( bmatched[0] == 2 && bmatched[1] == 3 ) mssmhbb.fill1DHistogram("jet_matching","matched_ranks",2.,mssmhbb.weight());
            }
            if ( ! mssmhbb.onlineBJetMatching({1,2,3},2)  )  continue;   // bjet trg matching (more inclusive)
            mssmhbb.fillJetHistograms("bjet_matching3");
            if ( ! mssmhbb.onlineBJetMatching({1,2},2)  )  continue;   // bjet trg matching (exclusive)
            mssmhbb.fillJetHistograms("bjet_matching2");
            mssmhbb.endSelection();
        } // end loop workflow 3
    }

} 

// end main
     



